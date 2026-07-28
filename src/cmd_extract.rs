//! # Extract subcommand
//!
//! The `extract` subcommand reads a FASTQ/A file, searches for subsequences
//! and writes matching records to a new FASTQ/A file. The output is written
//! to a new FASTQ/A file, with the file format determined by the input file.
//! Also, print detailed match information to stdout or a file if provided.
//!
//! The search algorithm is automatically selected based on the number of patterns
//! and their length. The BNDMq algorithm is used by default, but the user can
//! manually set the size of the _q_-grams. If the number of patterns is high or
//! the patterns are long, the Aho-Corasick algorithm is used.

use anyhow::{Context, Result};
use clap::{ArgAction, ArgGroup, Args, crate_name, crate_version};
use jiff::{Unit, Zoned};
use serde_json;

use std::collections::HashMap;
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::string::String;
use std::{env, fs};

use crossbeam_channel::{Receiver, Sender, bounded};
use paraseq::{Record, fastx};
use std::sync::Arc;

use crate::extract_processing::{
    ExtractSummary, FileSlot, IndexedResult, PipelineConfig, run_bounded_ordered_pipeline,
};
use crate::fastx_output::{FastxFormat, FastxRecordView, write_fastx_record};
use crate::helpers::{
    add_suffix_to_file_prefix, check_log_flag_conflict, error_if_directory,
    identify_uncompressed_type, parse_pattern_list, recommend_aho_corasick,
};
use crate::logger::{BufferedLogger, JsonLogger, append_json_log_fields, append_log_fields};
use crate::pattern_matching::PatternMatcher;

const DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE: usize = 1024;

#[derive(Debug, PartialEq, Eq)]
struct ThreadResolution {
    effective_total_threads: usize,
    auto_detected: bool,
    clamped: bool,
}

fn resolve_extract_threads_with_available(
    requested_threads: usize,
    available_threads: usize,
) -> ThreadResolution {
    let available_threads = available_threads.max(1);
    if requested_threads == 0 {
        return ThreadResolution {
            effective_total_threads: available_threads,
            auto_detected: true,
            clamped: false,
        };
    }

    ThreadResolution {
        effective_total_threads: requested_threads.min(available_threads),
        auto_detected: false,
        clamped: requested_threads > available_threads,
    }
}

fn resolve_extract_threads(requested_threads: usize) -> ThreadResolution {
    let available_threads = std::thread::available_parallelism()
        .map(|parallelism| parallelism.get())
        .unwrap_or(1);
    resolve_extract_threads_with_available(requested_threads, available_threads)
}

fn matching_thread_count(total_threads: usize) -> usize {
    total_threads.saturating_sub(1)
}

struct SingleRecordSetWorkChunk {
    start_index: u64,
    record_set: fastx::RecordSet,
    format: FastxFormat,
}

struct PairedRecordSetWorkChunk {
    start_index: u64,
    record_set_1: fastx::RecordSet,
    record_set_2: fastx::RecordSet,
    format_1: FastxFormat,
    format_2: FastxFormat,
}

struct SingleChunkResult {
    start_index: u64,
    record_count: usize,
    logs: ChunkLogs,
    summary: ExtractSummary,
    output: Vec<u8>,
}

impl IndexedResult for SingleChunkResult {
    fn index(&self) -> u64 {
        self.start_index
    }

    fn index_span(&self) -> u64 {
        self.record_count as u64
    }
}

struct PairedChunkResult {
    start_index: u64,
    record_count: usize,
    logs: ChunkLogs,
    summary: ExtractSummary,
    output_1: Vec<u8>,
    output_2: Vec<u8>,
}

impl IndexedResult for PairedChunkResult {
    fn index(&self) -> u64 {
        self.start_index
    }

    fn index_span(&self) -> u64 {
        self.record_count as u64
    }
}

fn record_set_len(record_set: &fastx::RecordSet) -> usize {
    match record_set {
        fastx::RecordSet::Fasta(records) => records.n_records(),
        fastx::RecordSet::Fastq(records) => records.n_records(),
    }
}

fn record_set_pool_size(config: PipelineConfig) -> usize {
    config.worker_count + config.work_queue_bound
}

#[derive(Default)]
struct ChunkLogs {
    plain: String,
    json: Vec<u8>,
    json_first: bool,
}

struct ChunkProcessor {
    matcher: Arc<PatternMatcher>,
    patterns: Arc<Vec<String>>,
    file_names: [String; 2],
    pattern_count: usize,
    logging_active: bool,
    plain_logging_active: bool,
    json_logging_active: bool,
    invert_match: bool,
    write_output: bool,
}

fn process_record_matches(
    processor: &ChunkProcessor,
    file_slot: FileSlot,
    record_id: &[u8],
    seq: &[u8],
    summary: &mut ExtractSummary,
    logs: &mut ChunkLogs,
) -> bool {
    if !processor.logging_active {
        return processor.matcher.find_any(seq);
    }

    summary.record_searched(seq.len());
    let mut matched = false;
    let file_name = match file_slot {
        FileSlot::SingleOrFirst => &processor.file_names[0],
        FileSlot::Second => &processor.file_names[1],
    };
    processor.matcher.for_each_match(seq, |hit| {
        matched = true;
        summary.pattern_hit(file_slot, hit.pattern_index);
        let pattern = &processor.patterns[hit.pattern_index];
        if processor.plain_logging_active {
            append_log_fields(&mut logs.plain, file_name, record_id, pattern, hit.position);
        }
        if processor.json_logging_active {
            append_json_log_fields(
                &mut logs.json,
                &mut logs.json_first,
                file_name,
                record_id,
                pattern,
                hit.position,
            );
        }
    });
    if matched {
        summary.record_hit(file_slot);
    }
    matched
}

fn process_borrowed_single_record_with_chunk_output<R: Record>(
    processor: &ChunkProcessor,
    record: R,
    format: FastxFormat,
    chunk: &mut SingleChunkResult,
) -> Result<()> {
    let seq = record.seq();
    let extracted = process_record_matches(
        processor,
        FileSlot::SingleOrFirst,
        record.id(),
        &seq,
        &mut chunk.summary,
        &mut chunk.logs,
    ) != processor.invert_match;
    if extracted {
        chunk.summary.extracted_records(1);
    }
    if processor.write_output && extracted {
        write_fastx_record(
            &mut chunk.output,
            FastxRecordView {
                id: record.id(),
                seq: &seq,
                qual: record.qual(),
                format,
            },
        )?;
    }
    Ok(())
}

fn process_borrowed_paired_record_with_chunk_output<R1: Record, R2: Record>(
    processor: &ChunkProcessor,
    record_1: R1,
    format_1: FastxFormat,
    record_2: R2,
    format_2: FastxFormat,
    chunk: &mut PairedChunkResult,
) -> Result<()> {
    let seq_1 = record_1.seq();
    let seq_2 = record_2.seq();
    let matched_1 = process_record_matches(
        processor,
        FileSlot::SingleOrFirst,
        record_1.id(),
        &seq_1,
        &mut chunk.summary,
        &mut chunk.logs,
    );
    let matched_2 = if !processor.logging_active && matched_1 {
        false
    } else {
        process_record_matches(
            processor,
            FileSlot::Second,
            record_2.id(),
            &seq_2,
            &mut chunk.summary,
            &mut chunk.logs,
        )
    };
    let extracted = (matched_1 || matched_2) != processor.invert_match;
    if extracted {
        chunk.summary.extracted_records(2);
    }
    if processor.write_output && extracted {
        write_fastx_record(
            &mut chunk.output_1,
            FastxRecordView {
                id: record_1.id(),
                seq: &seq_1,
                qual: record_1.qual(),
                format: format_1,
            },
        )?;
        write_fastx_record(
            &mut chunk.output_2,
            FastxRecordView {
                id: record_2.id(),
                seq: &seq_2,
                qual: record_2.qual(),
                format: format_2,
            },
        )?;
    }
    Ok(())
}

fn process_single_record_set(
    processor: &ChunkProcessor,
    start_index: u64,
    record_set: &fastx::RecordSet,
    format: FastxFormat,
) -> Result<SingleChunkResult> {
    let record_count = record_set_len(record_set);
    let mut chunk = SingleChunkResult {
        start_index,
        record_count,
        logs: ChunkLogs {
            json_first: true,
            ..ChunkLogs::default()
        },
        summary: ExtractSummary::new(processor.pattern_count),
        output: Vec::new(),
    };
    match record_set {
        fastx::RecordSet::Fasta(records) => {
            for record in records.iter() {
                let record = record.with_context(|| "Error during FASTQ/A record parsing.")?;
                process_borrowed_single_record_with_chunk_output(
                    processor, record, format, &mut chunk,
                )?;
            }
        }
        fastx::RecordSet::Fastq(records) => {
            for record in records.iter() {
                let record = record.with_context(|| "Error during FASTQ/A record parsing.")?;
                process_borrowed_single_record_with_chunk_output(
                    processor, record, format, &mut chunk,
                )?;
            }
        }
    }
    Ok(chunk)
}

fn process_pooled_single_record_set_chunk(
    processor: &ChunkProcessor,
    chunk: SingleRecordSetWorkChunk,
    pool_tx: &Sender<fastx::RecordSet>,
) -> Result<SingleChunkResult> {
    let SingleRecordSetWorkChunk {
        start_index,
        record_set,
        format,
    } = chunk;
    let result = process_single_record_set(processor, start_index, &record_set, format);
    pool_tx
        .send(record_set)
        .map_err(|_| anyhow::anyhow!("Single-end record-set pool closed during processing."))?;
    result
}

fn process_paired_record_iters<R1, R2, I1, I2>(
    processor: &ChunkProcessor,
    start_index: u64,
    capacity: usize,
    records: (I1, I2),
    formats: (FastxFormat, FastxFormat),
) -> Result<PairedChunkResult>
where
    R1: Record,
    R2: Record,
    I1: Iterator<Item = std::result::Result<R1, paraseq::Error>>,
    I2: Iterator<Item = std::result::Result<R2, paraseq::Error>>,
{
    let mut chunk = PairedChunkResult {
        start_index,
        record_count: capacity,
        logs: ChunkLogs {
            json_first: true,
            ..ChunkLogs::default()
        },
        summary: ExtractSummary::new(processor.pattern_count),
        output_1: Vec::new(),
        output_2: Vec::new(),
    };
    for (record_1, record_2) in records.0.zip(records.1) {
        let record_1 =
            record_1.with_context(|| "Error during FASTQ record parsing of first file.")?;
        let record_2 =
            record_2.with_context(|| "Error during FASTQ record parsing of second file.")?;
        process_borrowed_paired_record_with_chunk_output(
            processor, record_1, formats.0, record_2, formats.1, &mut chunk,
        )?;
    }
    Ok(chunk)
}

fn process_paired_record_sets(
    processor: &ChunkProcessor,
    start_index: u64,
    record_set_1: &fastx::RecordSet,
    record_set_2: &fastx::RecordSet,
    formats: (FastxFormat, FastxFormat),
) -> Result<PairedChunkResult> {
    let len_1 = record_set_len(record_set_1);
    let len_2 = record_set_len(record_set_2);
    if len_1 != len_2 {
        anyhow::bail!(
            "The two input files have a different number of records. Please provide valid paired-end read files."
        );
    }

    Ok(match (record_set_1, record_set_2) {
        (fastx::RecordSet::Fasta(records_1), fastx::RecordSet::Fasta(records_2)) => {
            process_paired_record_iters(
                processor,
                start_index,
                len_1,
                (records_1.iter(), records_2.iter()),
                formats,
            )?
        }
        (fastx::RecordSet::Fasta(records_1), fastx::RecordSet::Fastq(records_2)) => {
            process_paired_record_iters(
                processor,
                start_index,
                len_1,
                (records_1.iter(), records_2.iter()),
                formats,
            )?
        }
        (fastx::RecordSet::Fastq(records_1), fastx::RecordSet::Fasta(records_2)) => {
            process_paired_record_iters(
                processor,
                start_index,
                len_1,
                (records_1.iter(), records_2.iter()),
                formats,
            )?
        }
        (fastx::RecordSet::Fastq(records_1), fastx::RecordSet::Fastq(records_2)) => {
            process_paired_record_iters(
                processor,
                start_index,
                len_1,
                (records_1.iter(), records_2.iter()),
                formats,
            )?
        }
    })
}

fn process_pooled_paired_record_set_chunk(
    processor: &ChunkProcessor,
    chunk: PairedRecordSetWorkChunk,
    pool_tx: &Sender<(fastx::RecordSet, fastx::RecordSet)>,
) -> Result<PairedChunkResult> {
    let PairedRecordSetWorkChunk {
        start_index,
        record_set_1,
        record_set_2,
        format_1,
        format_2,
    } = chunk;
    let result = process_paired_record_sets(
        processor,
        start_index,
        &record_set_1,
        &record_set_2,
        (format_1, format_2),
    );
    pool_tx
        .send((record_set_1, record_set_2))
        .map_err(|_| anyhow::anyhow!("Paired-end record-set pool closed during processing."))?;
    result
}

fn consume_single_chunk_result(
    chunk: SingleChunkResult,
    writer: &mut Box<dyn io::Write>,
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
    summary: &mut ExtractSummary,
) -> Result<()> {
    logger.log_fragment(&chunk.logs.plain)?;
    if let Some(json_logger) = json_logger {
        json_logger.log_fragment(&chunk.logs.json)?;
    }
    summary.merge(&chunk.summary);
    if !chunk.output.is_empty() {
        io::Write::write_all(writer, &chunk.output)?;
    }
    Ok(())
}

fn consume_paired_chunk_result(
    chunk: PairedChunkResult,
    writers: (&mut Box<dyn io::Write>, &mut Box<dyn io::Write>),
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
    summary: &mut ExtractSummary,
) -> Result<()> {
    logger.log_fragment(&chunk.logs.plain)?;
    if let Some(json_logger) = json_logger {
        json_logger.log_fragment(&chunk.logs.json)?;
    }
    summary.merge(&chunk.summary);
    if !chunk.output_1.is_empty() {
        io::Write::write_all(writers.0, &chunk.output_1)?;
    }
    if !chunk.output_2.is_empty() {
        io::Write::write_all(writers.1, &chunk.output_2)?;
    }
    Ok(())
}

#[derive(Args)]
#[clap(group(
    ArgGroup::new("kmers")
        .required(true)
        .multiple(false)
        .args(&["kmer_seq", "kmer_file"]),
),
group(
    ArgGroup::new("algorithm")
        .required(false)
        .multiple(false)
        .args(&["q_size", "aho_corasick"]),
),
group(
    ArgGroup::new("logging")
        .required(false)
        .multiple(true)
        .args(&["out_log", "json_log"]),
),
group(
    ArgGroup::new("case-sensitivity")
        .required(false)
        .multiple(false)
        .args(&["case_insensitive", "lowercase", "uppercase"]),
),
group(
    ArgGroup::new("kmer-preprocessing")
        .required(false)
        .multiple(false)
        .args(&["canonical", "reverse_complement"])
))]
#[derive(Clone, Debug)]
pub struct CmdExtract {
    /// Input path for (compressed) FASTQ/A file.
    #[clap(short = 'i', long, short_alias = '1')]
    in_fastx: PathBuf,

    /// Input path for second FASTQ file (only for paired-end read processing).
    #[clap(short = '2', long, required = false)]
    in_fastq_2: Option<PathBuf>,

    /// Query sequences (accepts multiple sequences after the flag, separated by a space); if not provided, input path for file containing list of k-mers is required.
    #[clap(short = 's', long, num_args = 1..)]
    kmer_seq: Option<Vec<String>>,

    /// Input path for file containing list of k-mers, one per line (FASTA or plain text file; comment lines starting with '#'are ignored)
    #[clap(short = 'f', long)]
    kmer_file: Option<PathBuf>,

    /// Output file path for FASTQ/A file (extension derived from input file); if not provided, output is written to stdout.
    #[clap(short = 'o', long, required = false)]
    out_fastx: Option<PathBuf>,

    /// Also search for reverse complements of k-mers.
    #[clap(short = 'r', long, action(ArgAction::SetTrue), default_value("false"))]
    reverse_complement: bool,

    /// Search only for the canonical forms of k-mers.
    #[clap(short = 'c', long, action(ArgAction::SetTrue), default_value("false"))]
    canonical: bool,

    /// Print detailed match information to stdout if only the flag is passed, or to a file if a path is provided.
    #[clap(short = 'l', long, default_value(None), default_missing_value("STDOUT"), num_args = 0..=1, )]
    out_log: Option<PathBuf>,

    /// Write JSON log to stdout if only the flag is passed, or to a file if a path is provided.
    #[clap(short = 'j', long, default_value(None), default_missing_value("STDOUT"), num_args = 0..=1, )]
    json_log: Option<PathBuf>,

    /// Suppress output of found records (no records are written to a file or stdout); use if only matching statistics are of interest.
    #[clap(
        short = 'S',
        long,
        action(ArgAction::SetTrue),
        default_value("false"),
        conflicts_with("out_fastx"),
        requires("logging")
    )]
    suppress_output: bool,

    /// Invert the sense of matching, to select non-matching records.
    #[clap(short = 'v', long, action(ArgAction::SetTrue), default_value("false"))]
    invert_match: bool,

    /// Use case-insensitive matching. Always uses the Aho-Corasick algorithm.
    #[clap(short = 'I', long, action(ArgAction::SetTrue), default_value("false"))]
    case_insensitive: bool,

    /// Convert all input sequences to lowercase.
    #[clap(short = 'L', long, action(ArgAction::SetTrue), default_value("false"))]
    lowercase: bool,

    /// Convert all input sequences to uppercase.
    #[clap(short = 'U', long, action(ArgAction::SetTrue), default_value("false"))]
    uppercase: bool,

    /// Manually set size of q-grams to force the use of the BNDMq algorithm.
    #[clap(short = 'q', long, hide_short_help = true)]
    q_size: Option<usize>,

    /// Use Aho-Corasick to search for k-mers (best for lots of k-mers, or k-mers with more than 64 characters).
    #[clap(
        short = 'a',
        long,
        action(ArgAction::SetTrue),
        default_value("false"),
        hide_short_help = true
    )]
    aho_corasick: bool,

    /// Total number of processing threads. One thread reads input and the remaining threads match records. Use 0 to auto-detect available cores.
    #[clap(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Number of FASTA/Q records per parallel extract work chunk.
    #[clap(long, default_value_t = DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE, hide = true)]
    chunk_size: usize,
}

pub fn extract_records(args: CmdExtract) -> Result<()> {
    // Use helper for log flag conflict (not possible yet with `clap`)
    check_log_flag_conflict(
        &args.out_log,
        &args.json_log,
        &args.out_fastx,
        args.suppress_output,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    let mut args = args;

    let pattern_list = parse_pattern_list(
        &args.kmer_file,
        args.kmer_seq,
        args.reverse_complement,
        args.canonical,
        args.lowercase,
        args.uppercase,
    )
    .with_context(|| "Problem parsing pattern list.")?;

    // Case-insensitive matching always uses the Aho-Corasick algorithm
    if args.case_insensitive {
        args.aho_corasick = true;
    // Optimize search parameters only if user did not provide them
    } else if args.q_size.is_none() && !args.aho_corasick {
        args.aho_corasick = recommend_aho_corasick(&pattern_list)?;
    }

    // Set one of thre possible logging options:
    // 1) log to stdout,
    // 2) log to file,
    // 3) no logging
    let log_file = match &args.out_log {
        Some(path) => {
            if path
                .to_str()
                .ok_or_else(|| anyhow::anyhow!("Invalid log file path."))?
                == "STDOUT"
            {
                Some(Box::new(BufWriter::new(io::stdout())) as Box<dyn io::Write>)
            } else {
                let file = fs::File::create(path)
                    .with_context(|| format!("Problem creating log file: {}", path.display()))?;
                Some(Box::new(BufWriter::new(file)) as Box<dyn io::Write>)
            }
        }
        None => None,
    };

    //
    // -------------------- Initialization & Preprocessing --------------------
    //

    // Check if file paths point to directories and gets the file names
    error_if_directory(&args.in_fastx, "Record file path")?;
    let in_fastx_filename = args.in_fastx.file_name().unwrap().to_str().unwrap();
    let in_fastq_2_filename = match &args.in_fastq_2 {
        Some(p) => {
            error_if_directory(p, "Second read file path")?;
            p.file_name().unwrap().to_str().unwrap()
        }
        None => "",
    };

    // Activate logging if a log or JSON log file is provided
    let plain_logging_active = log_file.is_some();
    let json_logging_active = args.json_log.is_some();
    let logging_active = plain_logging_active || json_logging_active;

    // Initialize buffered logger with 8KB buffer
    let mut logger = BufferedLogger::new(log_file, 8192);

    // Initialize JSON logger if requested
    let mut json_logger = if let Some(json_path) = args.json_log.clone() {
        let writer: Box<dyn io::Write> = if json_path.to_str() == Some("STDOUT") {
            Box::new(io::stdout())
        } else {
            Box::new(fs::File::create(&json_path).with_context(|| {
                format!("Error creating JSON log file: {}", json_path.display())
            })?)
        };
        Some(JsonLogger::new(Some(writer), 8192))
    } else {
        None
    };

    // Log the list of patterns and header line
    if logging_active {
        // Write header section
        logger.write_header("#MerKurio extract log\n");
        logger.write_header(&format!("#{}\n", Zoned::now().round(Unit::Second)?));
        logger.write_header(&format!(
            "#Running {} version {}\n",
            crate_name!(),
            crate_version!()
        ));
        logger.write_header(&format!(
            "#Command line: {}\n",
            env::args().collect::<Vec<String>>().join(" ")
        ));
        logger.write_header(&format!(
            "#Searching for {} pattern{} {}\n",
            pattern_list.len(),
            if pattern_list.len() > 1 { "s" } else { "" },
            if args.invert_match {
                "(inverted matching)"
            } else {
                ""
            }
        ));
        logger.write_header("#\n#File\tRecord\tPattern\tPosition (zero-based)\n");
        logger.flush(); // Ensure header is written before records
    }

    let thread_resolution = resolve_extract_threads(args.threads);
    if thread_resolution.clamped {
        eprintln!(
            "Warning: requested {} extract threads, but only {} are available; using {} threads.",
            args.threads,
            thread_resolution.effective_total_threads,
            thread_resolution.effective_total_threads
        );
    }
    let total_threads = thread_resolution.effective_total_threads;
    let matching_threads = matching_thread_count(total_threads);
    if args.chunk_size == 0 {
        anyhow::bail!("Extract chunk size must be greater than zero.");
    }
    let parallel_chunk_size = args.chunk_size;

    let matcher = Arc::new(PatternMatcher::new(
        &pattern_list,
        args.aho_corasick,
        args.case_insensitive,
        args.q_size,
    )?);

    // Initialize counters for logging information
    let mut nb_records_tot = 0;
    let mut nb_bases = 0;
    let mut nb_hits_tot = (0, 0);
    let mut nb_records_hit = (0, 0);
    let mut nb_records_extracted = 0;
    let mut pattern_hit_counts: Vec<u32> = vec![0; pattern_list.len()];

    //
    // ------------------ Pattern Matching & Output Writing -------------------
    //

    //
    // If no second file is provided, process single file
    if args.in_fastq_2.is_none() {
        let mut reader = fastx::Reader::from_path(&args.in_fastx)
            .with_context(|| format!("Invalid FASTQ/A input path or file: {:?}", &args.in_fastx))?;
        let output_format = match reader.format() {
            fastx::Format::Fasta => FastxFormat::Fasta,
            fastx::Format::Fastq => FastxFormat::Fastq,
        };

        // Either write to file or stdout if no output path is provided;
        // the file format is determined by the input file
        let mut writer = match &args.out_fastx {
            Some(pathbuf) => {
                let mut out_path = pathbuf.clone();
                // Only add extension if not already present
                if out_path.extension().is_none() {
                    out_path = out_path
                        .with_extension(identify_uncompressed_type(&args.in_fastx).unwrap());
                }
                let path = Path::new(&out_path);
                let file = fs::File::create(path).with_context(|| {
                    format!("Error writing to output file; no such directory: {path:?}",)
                })?;
                Box::new(BufWriter::new(file)) as Box<dyn io::Write>
            }
            None => {
                let stdout = io::stdout();
                Box::new(BufWriter::new(stdout)) as Box<dyn io::Write>
            }
        };

        if total_threads > 1 {
            let write_output = !args.suppress_output;
            let chunk_processor = ChunkProcessor {
                matcher: Arc::clone(&matcher),
                patterns: Arc::new(pattern_list.clone()),
                file_names: [in_fastx_filename.to_string(), String::new()],
                pattern_count: pattern_list.len(),
                logging_active,
                plain_logging_active,
                json_logging_active,
                invert_match: args.invert_match,
                write_output,
            };
            let mut summary = ExtractSummary::new(pattern_list.len());
            let pipeline_config = PipelineConfig::new(matching_threads);
            let pool_size = record_set_pool_size(pipeline_config);
            let (record_pool_tx, record_pool_rx) = bounded::<fastx::RecordSet>(pool_size);
            for _ in 0..pool_size {
                record_pool_tx
                    .send(reader.new_record_set_with_size(parallel_chunk_size))
                    .map_err(|_| {
                        anyhow::anyhow!("Failed to initialize single-end record-set pool.")
                    })?;
            }
            let producer_record_pool_rx: Receiver<fastx::RecordSet> = record_pool_rx.clone();
            let worker_record_pool_tx = record_pool_tx.clone();
            run_bounded_ordered_pipeline(
                pipeline_config,
                move |work_tx| {
                    let mut start_index = 0_u64;
                    loop {
                        let mut record_set = producer_record_pool_rx.recv().map_err(|_| {
                            anyhow::anyhow!(
                                "Single-end record-set pool closed before parsing completed."
                            )
                        })?;
                        if !record_set
                            .fill(&mut reader)
                            .with_context(|| "Error during FASTQ/A record parsing.")?
                        {
                            break;
                        }
                        let len = record_set_len(&record_set);
                        let chunk = SingleRecordSetWorkChunk {
                            start_index,
                            record_set,
                            format: output_format,
                        };
                        work_tx.send(chunk).map_err(|_| {
                            anyhow::anyhow!(
                                "Pipeline work queue closed before all single-end work was sent."
                            )
                        })?;
                        start_index += len as u64;
                    }
                    Ok(())
                },
                move |chunk| {
                    process_pooled_single_record_set_chunk(
                        &chunk_processor,
                        chunk,
                        &worker_record_pool_tx,
                    )
                },
                |result| {
                    consume_single_chunk_result(
                        result,
                        &mut writer,
                        &mut logger,
                        &mut json_logger,
                        &mut summary,
                    )
                },
            )?;
            nb_records_tot = summary.nb_records_tot;
            nb_bases = summary.nb_bases;
            nb_hits_tot = summary.nb_hits_tot;
            nb_records_hit = summary.nb_records_hit;
            nb_records_extracted = summary.nb_records_extracted;
            pattern_hit_counts = summary.pattern_hit_counts;
        } else {
            // Iterate over FASTA/Q records and check for k-mer presence
            let mut record_set = reader.new_record_set();
            while record_set
                .fill(&mut reader)
                .with_context(|| "Error during FASTQ/A record parsing.")?
            {
                for record in record_set.iter() {
                    let record = record.with_context(|| "Error during FASTQ/A record parsing.")?;
                    let mut found_occ = false;
                    let seq = record.seq();

                    if logging_active {
                        nb_records_tot += 1;
                        nb_bases += seq.len();
                    }

                    // Report all matches when logging is active.
                    if logging_active {
                        matcher.for_each_match(&seq, |hit| {
                            found_occ = true;
                            logger.log_fields(
                                in_fastx_filename,
                                record.id(),
                                &pattern_list[hit.pattern_index],
                                hit.position,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastx_filename,
                                    record.id(),
                                    &pattern_list[hit.pattern_index],
                                    hit.position,
                                );
                            }
                            pattern_hit_counts[hit.pattern_index] += 1;
                            nb_hits_tot.0 += 1;
                        });
                        if found_occ {
                            nb_records_hit.0 += 1;
                        }
                    // Or use the fast first-match path when no per-match reporting is needed.
                    } else {
                        found_occ = matcher.find_any(&seq);
                    }

                    // Write record to file or stdout if any k-mer has been found
                    if found_occ != args.invert_match {
                        nb_records_extracted += 1;
                        if !args.suppress_output {
                            write_fastx_record(
                                &mut writer,
                                FastxRecordView {
                                    id: record.id(),
                                    seq: &seq,
                                    qual: record.qual(),
                                    format: output_format,
                                },
                            )?;
                        }
                    }
                }
            }
        }
    ////
    //// ---------------------- Handling Paired-End Reads ---------------------
    ////
    // If a second file is provided, process paired-end reads
    } else {
        let mut reader = fastx::Reader::from_path(&args.in_fastx)
            .with_context(|| format!("Invalid FASTQ/A input path or file: {:?}", &args.in_fastx))?;
        let output_format_1 = match reader.format() {
            fastx::Format::Fasta => FastxFormat::Fasta,
            fastx::Format::Fastq => FastxFormat::Fastq,
        };
        let mut reader_2 = fastx::Reader::from_path(args.in_fastq_2.clone().unwrap())
            .with_context(|| {
                format!(
                    "Invalid second FASTQ input path or file: {:?}",
                    &args.in_fastq_2
                )
            })?;
        let output_format_2 = match reader_2.format() {
            fastx::Format::Fasta => FastxFormat::Fasta,
            fastx::Format::Fastq => FastxFormat::Fastq,
        };

        // Either write to file or stdout if no output path is provided;
        // the file format is determined by the input file;
        // write to two files with _1 and _2 suffixes for paired-end reads
        let mut writer = match &args.out_fastx {
            Some(pathbuf) => {
                let pathbuf =
                    pathbuf.with_extension(identify_uncompressed_type(&args.in_fastx).unwrap());
                let pathbuf = add_suffix_to_file_prefix(&pathbuf, "_1");
                let path = Path::new(&pathbuf);
                let file = fs::File::create(path).with_context(|| {
                    format!("Error writing to paired-end file; no such directory: {path:?}",)
                })?;
                Box::new(BufWriter::new(file)) as Box<dyn io::Write>
            }
            None => {
                let stdout = io::stdout();
                Box::new(BufWriter::new(stdout)) as Box<dyn io::Write>
            }
        };
        let mut writer2 = match args.out_fastx {
            Some(pathbuf) => {
                let pathbuf =
                    pathbuf.with_extension(identify_uncompressed_type(&args.in_fastx).unwrap());
                let pathbuf = add_suffix_to_file_prefix(&pathbuf, "_2");
                let path = Path::new(&pathbuf);
                let file = fs::File::create(path).with_context(|| {
                    format!("Error writing second paired-end file; no such directory: {path:?}",)
                })?;
                Box::new(BufWriter::new(file)) as Box<dyn io::Write>
            }
            None => {
                let stdout = io::stdout();
                Box::new(BufWriter::new(stdout)) as Box<dyn io::Write>
            }
        };

        if total_threads > 1 {
            let write_output = !args.suppress_output;
            let chunk_processor = ChunkProcessor {
                matcher: Arc::clone(&matcher),
                patterns: Arc::new(pattern_list.clone()),
                file_names: [
                    in_fastx_filename.to_string(),
                    in_fastq_2_filename.to_string(),
                ],
                pattern_count: pattern_list.len(),
                logging_active,
                plain_logging_active,
                json_logging_active,
                invert_match: args.invert_match,
                write_output,
            };
            let mut summary = ExtractSummary::new(pattern_list.len());
            let pipeline_config = PipelineConfig::new(matching_threads);
            let pool_size = record_set_pool_size(pipeline_config);
            let (record_pool_tx, record_pool_rx) =
                bounded::<(fastx::RecordSet, fastx::RecordSet)>(pool_size);
            for _ in 0..pool_size {
                record_pool_tx
                    .send((
                        reader.new_record_set_with_size(parallel_chunk_size),
                        reader_2.new_record_set_with_size(parallel_chunk_size),
                    ))
                    .map_err(|_| {
                        anyhow::anyhow!("Failed to initialize paired-end record-set pool.")
                    })?;
            }
            let producer_record_pool_rx: Receiver<(fastx::RecordSet, fastx::RecordSet)> =
                record_pool_rx.clone();
            let worker_record_pool_tx = record_pool_tx.clone();
            run_bounded_ordered_pipeline(
                pipeline_config,
                move |work_tx| {
                    let mut start_index = 0_u64;

                    loop {
                        let (mut record_set_1, mut record_set_2) =
                            producer_record_pool_rx.recv().map_err(|_| {
                                anyhow::anyhow!(
                                    "Paired-end record-set pool closed before parsing completed."
                                )
                            })?;
                        let filled_1 = record_set_1
                            .fill(&mut reader)
                            .with_context(|| "Error during FASTQ record parsing of first file.")?;
                        let filled_2 = record_set_2
                            .fill(&mut reader_2)
                            .with_context(|| "Error during FASTQ record parsing of second file.")?;

                        match (filled_1, filled_2) {
                            (false, false) => break,
                            (true, true) => {}
                            _ => {
                                anyhow::bail!(
                                    "The two input files have a different number of records. Please provide valid paired-end read files."
                                );
                            }
                        }

                        let len_1 = record_set_len(&record_set_1);
                        let len_2 = record_set_len(&record_set_2);
                        if len_1 != len_2 {
                            anyhow::bail!(
                                "The two input files have a different number of records. Please provide valid paired-end read files."
                            );
                        }

                        let chunk = PairedRecordSetWorkChunk {
                            start_index,
                            record_set_1,
                            record_set_2,
                            format_1: output_format_1,
                            format_2: output_format_2,
                        };
                        work_tx.send(chunk).map_err(|_| {
                            anyhow::anyhow!(
                                "Pipeline work queue closed before all paired-end work was sent."
                            )
                        })?;
                        start_index += len_1 as u64;
                    }
                    Ok(())
                },
                move |chunk| {
                    process_pooled_paired_record_set_chunk(
                        &chunk_processor,
                        chunk,
                        &worker_record_pool_tx,
                    )
                },
                |result| {
                    consume_paired_chunk_result(
                        result,
                        (&mut writer, &mut writer2),
                        &mut logger,
                        &mut json_logger,
                        &mut summary,
                    )
                },
            )?;
            nb_records_tot = summary.nb_records_tot;
            nb_bases = summary.nb_bases;
            nb_hits_tot = summary.nb_hits_tot;
            nb_records_hit = summary.nb_records_hit;
            nb_records_extracted = summary.nb_records_extracted;
            pattern_hit_counts = summary.pattern_hit_counts;
        } else {
            // Iterate over FASTQ records and check for k-mer presence
            let mut record_set_1 = reader.new_record_set();
            let mut record_set_2 = reader_2.new_record_set();
            loop {
                let filled_1 = record_set_1
                    .fill(&mut reader)
                    .with_context(|| "Error during FASTQ record parsing of first file.")?;
                let filled_2 = record_set_2
                    .fill(&mut reader_2)
                    .with_context(|| "Error during FASTQ record parsing of second file.")?;

                match (filled_1, filled_2) {
                    (false, false) => break,
                    (true, true) => {}
                    _ => {
                        anyhow::bail!(
                            "The two input files have a different number of records. Please provide valid paired-end read files."
                        );
                    }
                }

                if record_set_len(&record_set_1) != record_set_len(&record_set_2) {
                    anyhow::bail!(
                        "The two input files have a different number of records. Please provide valid paired-end read files."
                    );
                }

                for (record_1, record_2) in record_set_1.iter().zip(record_set_2.iter()) {
                    let record_1 = record_1
                        .with_context(|| "Error during FASTQ record parsing of first file.")?;
                    let record_2 = record_2
                        .with_context(|| "Error during FASTQ record parsing of second file.")?;
                    let mut found_occ = false;
                    let seq_1 = record_1.seq();
                    let seq_2 = record_2.seq();

                    if logging_active {
                        nb_records_tot += 2;
                        nb_bases += seq_1.len();
                        nb_bases += seq_2.len();
                    }

                    // Report all matches when logging is active.
                    if logging_active {
                        let mut record_hit: (usize, usize) = (0, 0);
                        matcher.for_each_match(&seq_1, |hit| {
                            logger.log_fields(
                                in_fastx_filename,
                                record_1.id(),
                                &pattern_list[hit.pattern_index],
                                hit.position,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastx_filename,
                                    record_1.id(),
                                    &pattern_list[hit.pattern_index],
                                    hit.position,
                                );
                            }
                            pattern_hit_counts[hit.pattern_index] += 1;
                            record_hit.0 = 1;
                            nb_hits_tot.0 += 1;
                            found_occ = true;
                        });
                        matcher.for_each_match(&seq_2, |hit| {
                            logger.log_fields(
                                in_fastq_2_filename,
                                record_2.id(),
                                &pattern_list[hit.pattern_index],
                                hit.position,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastq_2_filename,
                                    record_2.id(),
                                    &pattern_list[hit.pattern_index],
                                    hit.position,
                                );
                            }
                            pattern_hit_counts[hit.pattern_index] += 1;
                            record_hit.1 = 1;
                            nb_hits_tot.1 += 1;
                            found_occ = true;
                        });
                        nb_records_hit.0 += record_hit.0;
                        nb_records_hit.1 += record_hit.1;
                    // Or use the fast first-match path when no per-match reporting is needed.
                    } else {
                        found_occ = matcher.find_any(&seq_1) || matcher.find_any(&seq_2);
                    }

                    // Write records to file or stdout if any patterns have been matched
                    if found_occ != args.invert_match {
                        nb_records_extracted += 2;
                        if !args.suppress_output {
                            write_fastx_record(
                                &mut writer,
                                FastxRecordView {
                                    id: record_1.id(),
                                    seq: &seq_1,
                                    qual: record_1.qual(),
                                    format: output_format_1,
                                },
                            )?;
                            write_fastx_record(
                                &mut writer2,
                                FastxRecordView {
                                    id: record_2.id(),
                                    seq: &seq_2,
                                    qual: record_2.qual(),
                                    format: output_format_2,
                                },
                            )?;
                        }
                    }
                }
            }
        }
    }

    // Log summary statistics as plain text and/or JSON
    if logging_active {
        logger.flush();
        let nb_patterns_found = pattern_hit_counts
            .iter()
            .filter(|&&count| count > 0)
            .count();
        let nb_patterns_found_percentage =
            nb_patterns_found as f64 / pattern_hit_counts.len() as f64 * 100.0;
        logger.write_header(&format!(
            "#\n#Number of patterns found: {}/{} ({:.2} %)\n",
            nb_patterns_found,
            pattern_hit_counts.len(),
            nb_patterns_found_percentage,
        ));
        logger.write_header("#Pattern\tCount\n");
        for (pattern, count) in pattern_list.iter().zip(pattern_hit_counts.iter()) {
            logger.write_header(&format!("#{pattern}\t{count}\n"));
        }
        logger.write_header(&format!(
            "#\n#Total number of records searched: {nb_records_tot}\n"
        ));
        logger.write_header(&format!(
            "#Total number of characters searched: {nb_bases}\n"
        ));
        logger.write_header(&format!(
            "#Total number of hits: {}\n",
            nb_hits_tot.0 + nb_hits_tot.1
        ));
        logger.write_header(&format!(
            "#Number of distinct records with a hit: {}\n",
            nb_records_hit.0 + nb_records_hit.1
        ));
        if args.in_fastq_2.is_some() {
            logger.write_header(&format!(
                "#\n#Total number of hits in file 1: {}\n",
                nb_hits_tot.0
            ));
            logger.write_header(&format!(
                "#Total number of hits in file 2: {}\n",
                nb_hits_tot.1
            ));
            logger.write_header(&format!(
                "#Number of distinct records with a hit in file 1: {}\n",
                nb_records_hit.0
            ));
            logger.write_header(&format!(
                "#Number of distinct records with a hit in file 2: {}\n",
                nb_records_hit.1
            ));
            logger.write_header(&format!(
                "#Total number of extracted records: {nb_records_extracted}\n"
            ));
        }
        logger.flush();
    }

    // Finalize JSON log if active
    if let Some(jl) = json_logger {
        let input_files_json = serde_json::json!({
            "kmer_file": args.kmer_file.as_ref().map(|p| p.to_string_lossy().to_string()),
            "record_file_1": in_fastx_filename,
            "record_file_2": if args.in_fastq_2.is_some() { Some(in_fastq_2_filename) } else { None },
        });
        let pattern_hit_counts_map: HashMap<String, u32> = pattern_list
            .iter()
            .cloned()
            .zip(pattern_hit_counts.iter().copied())
            .collect();
        let meta_information = serde_json::json!({
            "program": crate_name!(),
            "version": crate_version!(),
            "timestamp": Zoned::now().round(Unit::Second).unwrap(),
            "subcommand": "extract",
            "command_line": env::args().collect::<Vec<String>>(),
            "search_algorithm": matcher.algorithm_name(),
            "inverted_matching": args.invert_match,
            "case_insensitive": args.case_insensitive,
            "input_files": input_files_json,
        });
        let summary_statistics = serde_json::json!({
            "number_of_patterns_searched": pattern_list.len(),
            "number_of_patterns_found": pattern_hit_counts.iter().filter(|&&count| count > 0).count(),
            "number_of_records_searched": nb_records_tot,
            "number_of_characters_searched": nb_bases,
            "number_of_matches": nb_hits_tot.0 + nb_hits_tot.1,
            "number_of_distinct_records_with_a_hit": nb_records_hit.0 + nb_records_hit.1,
        });
        let paired_end_stats = serde_json::json!({
            "searching_paired_end_reads": args.in_fastq_2.is_some(),
            "number_of_hits_in_file_1": nb_hits_tot.0,
            "number_of_hits_in_file_2": if args.in_fastq_2.is_some() { Some(nb_hits_tot.1) } else { None },
            "number_of_distinct_records_with_a_hit_in_file_1": nb_records_hit.0,
            "number_of_distinct_records_with_a_hit_in_file_2": if args.in_fastq_2.is_some() { Some(nb_records_hit.1) } else { None },
            "number_of_extracted_records": nb_records_extracted,
        });
        jl.finalize(
            &meta_information,
            &serde_json::json!(pattern_hit_counts_map),
            &summary_statistics,
            Some(&paired_end_stats),
        );
    }

    Ok(())
}

//
// ---------------------------------- Tests ----------------------------------
//

#[cfg(test)]
#[path = "cmd_extract/tests.rs"]
mod tests;
