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

use paraseq::{Record, fastx};
use std::sync::Arc;

use crate::extract_processing::{
    ExtractProcessor, ExtractSummary, FileSlot, OutputRecord, PairedResult, PipelineConfig,
    RecordHit, RecordInput, SingleResult, run_bounded_ordered_pipeline_with_producer,
};
use crate::fastx_output::{FastxFormat, FastxRecordView, write_fastx_record};
use crate::helpers::{
    add_suffix_to_file_prefix, check_log_flag_conflict, error_if_directory,
    identify_uncompressed_type, parse_pattern_list, recommend_aho_corasick,
};
use crate::logger::{BufferedLogger, JsonLogger};
use crate::pattern_matching::PatternMatcher;

const EXTRACT_PARALLEL_CHUNK_SIZE: usize = 512;

#[derive(Debug, PartialEq, Eq)]
struct ThreadResolution {
    effective_threads: usize,
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
            effective_threads: available_threads,
            auto_detected: true,
            clamped: false,
        };
    }

    ThreadResolution {
        effective_threads: requested_threads.min(available_threads),
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

#[derive(Debug)]
struct OwnedRecordInput {
    index: u64,
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
    format: FastxFormat,
}

#[derive(Debug)]
struct SingleWorkChunk {
    records: Vec<OwnedRecordInput>,
}

#[derive(Debug)]
struct OwnedPairedInput {
    index: u64,
    read_1: OwnedRecordInput,
    read_2: OwnedRecordInput,
}

#[derive(Debug)]
struct PairedWorkChunk {
    pairs: Vec<OwnedPairedInput>,
}

fn process_single_work_chunk(
    processor: &ExtractProcessor,
    chunk: SingleWorkChunk,
) -> Vec<SingleResult> {
    chunk
        .records
        .iter()
        .map(|record| {
            processor.process_single_record(
                record.index,
                RecordInput {
                    id: &record.id,
                    seq: &record.seq,
                    qual: record.qual.as_deref(),
                    format: record.format,
                },
            )
        })
        .collect()
}

fn process_paired_work_chunk(
    processor: &ExtractProcessor,
    chunk: PairedWorkChunk,
) -> Vec<PairedResult> {
    chunk
        .pairs
        .iter()
        .map(|pair| {
            processor.process_paired_record(
                pair.index,
                RecordInput {
                    id: &pair.read_1.id,
                    seq: &pair.read_1.seq,
                    qual: pair.read_1.qual.as_deref(),
                    format: pair.read_1.format,
                },
                RecordInput {
                    id: &pair.read_2.id,
                    seq: &pair.read_2.seq,
                    qual: pair.read_2.qual.as_deref(),
                    format: pair.read_2.format,
                },
            )
        })
        .collect()
}

fn output_record_view(record: &OutputRecord) -> FastxRecordView<'_> {
    FastxRecordView {
        id: &record.id,
        seq: &record.seq,
        qual: record.qual.as_deref(),
        format: record.format,
    }
}

fn log_record_hit(
    hit: &RecordHit,
    file_1_name: &str,
    file_2_name: &str,
    pattern_list: &[String],
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
) {
    let file_name = match hit.file_slot {
        FileSlot::SingleOrFirst => file_1_name,
        FileSlot::Second => file_2_name,
    };
    logger.log_fields(
        file_name,
        &hit.record_id,
        &pattern_list[hit.pattern_index],
        hit.position,
    );
    if let Some(jl) = json_logger {
        jl.log_fields(
            file_name,
            &hit.record_id,
            &pattern_list[hit.pattern_index],
            hit.position,
        );
    }
}

fn consume_single_result(
    result: SingleResult,
    file_name: &str,
    pattern_list: &[String],
    writer: &mut Box<dyn io::Write>,
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
    summary: &mut ExtractSummary,
) -> Result<()> {
    for hit in &result.hits {
        log_record_hit(hit, file_name, "", pattern_list, logger, json_logger);
    }
    if let Some(output) = &result.output {
        write_fastx_record(writer, output_record_view(output))?;
    }
    summary.merge_single(&result);
    Ok(())
}

fn consume_paired_result(
    result: PairedResult,
    file_1_name: &str,
    file_2_name: &str,
    pattern_list: &[String],
    writer_1: &mut Box<dyn io::Write>,
    writer_2: &mut Box<dyn io::Write>,
    logger: &mut BufferedLogger,
    json_logger: &mut Option<JsonLogger>,
    summary: &mut ExtractSummary,
) -> Result<()> {
    for hit in &result.hits {
        log_record_hit(
            hit,
            file_1_name,
            file_2_name,
            pattern_list,
            logger,
            json_logger,
        );
    }
    if let Some(output) = &result.output_1 {
        write_fastx_record(writer_1, output_record_view(output))?;
    }
    if let Some(output) = &result.output_2 {
        write_fastx_record(writer_2, output_record_view(output))?;
    }
    summary.merge_paired(&result);
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

    /// Number of worker threads for extract pattern matching. Use 0 to auto-detect available cores.
    #[clap(short = 't', long, default_value_t = 1)]
    threads: usize,
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
    let logging_active = log_file.is_some() || args.json_log.is_some();

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
            "Warning: requested {} extract worker threads, but only {} are available; using {} threads.",
            args.threads, thread_resolution.effective_threads, thread_resolution.effective_threads
        );
    }
    let worker_threads = thread_resolution.effective_threads;

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

        if worker_threads > 1 {
            let processor = ExtractProcessor::new(
                Arc::clone(&matcher),
                pattern_list.len(),
                logging_active,
                args.suppress_output,
                args.invert_match,
            );
            let worker_processor = processor.clone();
            let mut summary = ExtractSummary::new(pattern_list.len());
            run_bounded_ordered_pipeline_with_producer(
                PipelineConfig::new(worker_threads),
                move |work_tx| {
                    let mut record_set = reader.new_record_set();
                    let mut chunk = SingleWorkChunk {
                        records: Vec::with_capacity(EXTRACT_PARALLEL_CHUNK_SIZE),
                    };
                    let mut record_index = 0_u64;
                    while record_set
                        .fill(&mut reader)
                        .with_context(|| "Error during FASTQ/A record parsing.")?
                    {
                        for record in record_set.iter() {
                            let record =
                                record.with_context(|| "Error during FASTQ/A record parsing.")?;
                            let seq = record.seq();
                            chunk.records.push(OwnedRecordInput {
                                index: record_index,
                                id: record.id().to_vec(),
                                seq: seq.into_owned(),
                                qual: record.qual().map(|qual| qual.to_vec()),
                                format: output_format,
                            });
                            record_index += 1;
                            if chunk.records.len() == EXTRACT_PARALLEL_CHUNK_SIZE {
                                work_tx.send(chunk).map_err(|_| {
                                    anyhow::anyhow!(
                                        "Pipeline work queue closed before all single-end work was sent."
                                    )
                                })?;
                                chunk = SingleWorkChunk {
                                    records: Vec::with_capacity(EXTRACT_PARALLEL_CHUNK_SIZE),
                                };
                            }
                        }
                    }
                    if !chunk.records.is_empty() {
                        work_tx.send(chunk).map_err(|_| {
                            anyhow::anyhow!(
                                "Pipeline work queue closed before all single-end work was sent."
                            )
                        })?;
                    }
                    Ok(())
                },
                move |chunk| Ok(process_single_work_chunk(&worker_processor, chunk)),
                |result| {
                    consume_single_result(
                        result,
                        in_fastx_filename,
                        &pattern_list,
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

        if worker_threads > 1 {
            let processor = ExtractProcessor::new(
                Arc::clone(&matcher),
                pattern_list.len(),
                logging_active,
                args.suppress_output,
                args.invert_match,
            );
            let worker_processor = processor.clone();
            let mut summary = ExtractSummary::new(pattern_list.len());
            run_bounded_ordered_pipeline_with_producer(
                PipelineConfig::new(worker_threads),
                move |work_tx| {
                    let mut record_set_1 = reader.new_record_set();
                    let mut record_set_2 = reader_2.new_record_set();
                    let mut chunk = PairedWorkChunk {
                        pairs: Vec::with_capacity(EXTRACT_PARALLEL_CHUNK_SIZE),
                    };
                    let mut pair_index = 0_u64;

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

                        let records_1 = record_set_1
                            .iter()
                            .collect::<std::result::Result<Vec<_>, _>>()
                            .with_context(|| "Error during FASTQ record parsing of first file.")?;
                        let records_2 = record_set_2
                            .iter()
                            .collect::<std::result::Result<Vec<_>, _>>()
                            .with_context(|| "Error during FASTQ record parsing of second file.")?;

                        if records_1.len() != records_2.len() {
                            anyhow::bail!(
                                "The two input files have a different number of records. Please provide valid paired-end read files."
                            );
                        }

                        for (record_1, record_2) in records_1.into_iter().zip(records_2.into_iter())
                        {
                            let seq_1 = record_1.seq();
                            let seq_2 = record_2.seq();
                            chunk.pairs.push(OwnedPairedInput {
                                index: pair_index,
                                read_1: OwnedRecordInput {
                                    index: pair_index,
                                    id: record_1.id().to_vec(),
                                    seq: seq_1.into_owned(),
                                    qual: record_1.qual().map(|qual| qual.to_vec()),
                                    format: output_format_1,
                                },
                                read_2: OwnedRecordInput {
                                    index: pair_index,
                                    id: record_2.id().to_vec(),
                                    seq: seq_2.into_owned(),
                                    qual: record_2.qual().map(|qual| qual.to_vec()),
                                    format: output_format_2,
                                },
                            });
                            pair_index += 1;
                            if chunk.pairs.len() == EXTRACT_PARALLEL_CHUNK_SIZE {
                                work_tx.send(chunk).map_err(|_| {
                                    anyhow::anyhow!(
                                        "Pipeline work queue closed before all paired-end work was sent."
                                    )
                                })?;
                                chunk = PairedWorkChunk {
                                    pairs: Vec::with_capacity(EXTRACT_PARALLEL_CHUNK_SIZE),
                                };
                            }
                        }
                    }

                    if !chunk.pairs.is_empty() {
                        work_tx.send(chunk).map_err(|_| {
                            anyhow::anyhow!(
                                "Pipeline work queue closed before all paired-end work was sent."
                            )
                        })?;
                    }
                    Ok(())
                },
                move |chunk| Ok(process_paired_work_chunk(&worker_processor, chunk)),
                |result| {
                    consume_paired_result(
                        result,
                        in_fastx_filename,
                        in_fastq_2_filename,
                        &pattern_list,
                        &mut writer,
                        &mut writer2,
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

                let records_1 = record_set_1
                    .iter()
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .with_context(|| "Error during FASTQ record parsing of first file.")?;
                let records_2 = record_set_2
                    .iter()
                    .collect::<std::result::Result<Vec<_>, _>>()
                    .with_context(|| "Error during FASTQ record parsing of second file.")?;

                if records_1.len() != records_2.len() {
                    anyhow::bail!(
                        "The two input files have a different number of records. Please provide valid paired-end read files."
                    );
                }

                for (record_1, record_2) in records_1.into_iter().zip(records_2.into_iter()) {
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
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    struct SingleExtractRun {
        _temp_dir: tempfile::TempDir,
        out_fastx: PathBuf,
        out_log: PathBuf,
        out_json: PathBuf,
    }

    #[derive(Clone)]
    struct SingleExtractOptions {
        input: PathBuf,
        patterns: Vec<String>,
        threads: usize,
        aho_corasick: bool,
        case_insensitive: bool,
        invert_match: bool,
        suppress_output: bool,
        reverse_complement: bool,
        log: bool,
        json: bool,
    }

    impl SingleExtractOptions {
        fn simple(patterns: Vec<String>, threads: usize) -> Self {
            Self {
                input: PathBuf::from("tests/fixtures/input/simple.fasta"),
                patterns,
                threads,
                aho_corasick: false,
                case_insensitive: false,
                invert_match: false,
                suppress_output: false,
                reverse_complement: true,
                log: true,
                json: true,
            }
        }
    }

    fn run_single_extract(options: SingleExtractOptions) -> Result<SingleExtractRun> {
        let temp_dir = tempfile::tempdir()?;
        let out_fastx = temp_dir.path().join("out.fasta");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        let args = CmdExtract {
            in_fastx: options.input,
            in_fastq_2: None,
            kmer_seq: Some(options.patterns),
            kmer_file: None,
            out_fastx: if options.suppress_output {
                None
            } else {
                Some(out_fastx.clone())
            },
            q_size: None,
            aho_corasick: options.aho_corasick,
            reverse_complement: options.reverse_complement,
            canonical: false,
            out_log: options.log.then_some(out_log.clone()),
            suppress_output: options.suppress_output,
            json_log: options.json.then_some(out_json.clone()),
            invert_match: options.invert_match,
            case_insensitive: options.case_insensitive,
            lowercase: false,
            uppercase: false,
            threads: options.threads,
        };

        extract_records(args)?;

        Ok(SingleExtractRun {
            _temp_dir: temp_dir,
            out_fastx,
            out_log,
            out_json,
        })
    }

    /// Compare FASTA output with fixture
    fn compare_fasta_output(actual_path: &Path, expected_path: &str) -> Result<()> {
        let expected = fs::read_to_string(expected_path)?;
        let actual = fs::read_to_string(actual_path)?;
        assert_eq!(expected, actual, "FASTA output does not match fixture");
        Ok(())
    }

    fn compare_text_files(actual_path: &Path, expected_path: &Path) -> Result<()> {
        let expected = fs::read_to_string(expected_path)?;
        let actual = fs::read_to_string(actual_path)?;
        assert_eq!(expected, actual);
        Ok(())
    }

    fn compare_log_outputs(actual_path: &Path, expected_path: &Path) -> Result<()> {
        let expected_log = fs::read_to_string(expected_path)?;
        let actual_log = fs::read_to_string(actual_path)?;
        let expected_lines: Vec<&str> = expected_log.lines().skip(4).collect();
        let actual_lines: Vec<&str> = actual_log.lines().skip(4).collect();
        assert_eq!(expected_lines, actual_lines);
        Ok(())
    }

    fn compare_json_outputs(actual_path: &Path, expected_path: &Path) -> Result<()> {
        let expected_json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(expected_path)?)?;
        let actual_json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(actual_path)?)?;

        assert_eq!(
            expected_json["matching_records"],
            actual_json["matching_records"]
        );
        assert_eq!(
            expected_json["summary_statistics"],
            actual_json["summary_statistics"]
        );
        assert_eq!(
            expected_json["paired_end_reads_statistics"],
            actual_json["paired_end_reads_statistics"]
        );
        assert_eq!(
            expected_json["pattern_hit_counts"],
            actual_json["pattern_hit_counts"]
        );
        assert_eq!(
            expected_json["meta_information"]["search_algorithm"],
            actual_json["meta_information"]["search_algorithm"]
        );
        assert_eq!(
            expected_json["meta_information"]["inverted_matching"],
            actual_json["meta_information"]["inverted_matching"]
        );
        assert_eq!(
            expected_json["meta_information"]["case_insensitive"],
            actual_json["meta_information"]["case_insensitive"]
        );
        Ok(())
    }

    fn assert_single_runs_match(
        actual: &SingleExtractRun,
        expected: &SingleExtractRun,
    ) -> Result<()> {
        compare_text_files(&actual.out_fastx, &expected.out_fastx)?;
        compare_log_outputs(&actual.out_log, &expected.out_log)?;
        compare_json_outputs(&actual.out_json, &expected.out_json)?;
        Ok(())
    }

    /// Compare log output with fixture, ignoring metadata
    fn compare_log_output(actual_path: &Path, expected_path: &str) -> Result<()> {
        let expected_log = fs::read_to_string(expected_path)?;
        let actual_log = fs::read_to_string(actual_path)?;

        // Split into lines and compare, skipping metadata lines
        let expected_lines: Vec<&str> = expected_log.lines().collect();
        let actual_lines: Vec<&str> = actual_log.lines().collect();

        // Skip first 4 lines (metadata) and compare the rest
        let expected_content = &expected_lines[4..];
        let actual_content = &actual_lines[4..];

        // Find section boundaries in expected content
        let mut section_boundaries = Vec::new();
        let mut in_match_section = false;
        let mut in_pattern_section = false;

        for (i, line) in expected_content.iter().enumerate() {
            if line.starts_with('#') {
                if line.contains("Pattern\tCount") {
                    in_pattern_section = true;
                    section_boundaries.push(i);
                } else if line.contains("Number of patterns found") {
                    in_match_section = false;
                    section_boundaries.push(i);
                } else if !in_match_section && !in_pattern_section {
                    in_match_section = true;
                    section_boundaries.push(i);
                }
            }
        }
        section_boundaries.push(expected_content.len());

        // Compare pattern count and header separator
        assert_eq!(
            expected_content[0], actual_content[0],
            "Log pattern count mismatch"
        );
        assert_eq!(
            expected_content[1], actual_content[1],
            "Log header separator mismatch"
        );

        // Compare column headers
        assert_eq!(
            expected_content[2], actual_content[2],
            "Log column header mismatch"
        );

        // Compare match records (between header and pattern count summary)
        let match_start = 3;
        let match_end = section_boundaries[1];
        for i in match_start..match_end {
            assert_eq!(
                expected_content[i],
                actual_content[i],
                "Log match record mismatch at line {}",
                i + 5
            );
        }

        // Compare pattern count summary header and separator
        assert_eq!(
            expected_content[match_end], actual_content[match_end],
            "Log pattern count summary header mismatch"
        );
        assert_eq!(
            expected_content[match_end + 1],
            actual_content[match_end + 1],
            "Log pattern count header mismatch"
        );

        // Compare pattern counts
        let pattern_start = match_end + 2;
        let pattern_end = section_boundaries[2];
        for i in pattern_start..pattern_end {
            assert_eq!(
                expected_content[i],
                actual_content[i],
                "Log pattern count mismatch at line {}",
                i + 5
            );
        }

        // Compare summary statistics
        let stats_start = pattern_end + 1;
        let stats_end = section_boundaries[3];
        for i in stats_start..stats_end {
            assert_eq!(
                expected_content[i],
                actual_content[i],
                "Log summary statistic mismatch at line {}",
                i + 5
            );
        }

        Ok(())
    }

    /// Compare JSON output with fixture, mostly ignoring metadata
    fn compare_json_output(actual_path: &Path, expected_path: &str) -> Result<()> {
        let expected_json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(expected_path)?)?;
        let actual_json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(actual_path)?)?;

        // Compare all fields except meta_information
        assert_eq!(
            expected_json["matching_records"], actual_json["matching_records"],
            "JSON matching records mismatch"
        );
        assert_eq!(
            expected_json["summary_statistics"], actual_json["summary_statistics"],
            "JSON summary statistics mismatch"
        );
        assert_eq!(
            expected_json["paired_end_reads_statistics"],
            actual_json["paired_end_reads_statistics"],
            "JSON paired-end reads statistics mismatch"
        );
        assert_eq!(
            expected_json["pattern_hit_counts"], actual_json["pattern_hit_counts"],
            "JSON pattern hit counts mismatch"
        );

        // Compare specific meta_information fields
        assert_eq!(
            expected_json["meta_information"]["search_algorithm"],
            actual_json["meta_information"]["search_algorithm"],
            "JSON search algorithm mismatch"
        );
        assert_eq!(
            expected_json["meta_information"]["inverted_matching"],
            actual_json["meta_information"]["inverted_matching"],
            "JSON inverted matching mismatch"
        );
        assert_eq!(
            expected_json["meta_information"]["case_insensitive"],
            actual_json["meta_information"]["case_insensitive"],
            "JSON case insensitive mismatch"
        );

        Ok(())
    }

    // Compare with simple nucleotide FASTA file, including the reverse complement
    // Corresponds to: cargo run -- extract -i tests/fixtures/input/simple.fasta -r -s ACG -o tests/fixtures/extract/simple.extracted.fasta -l tests/fixtures/extract/simple.log -j tests/fixtures/extract/simple.json
    #[test]
    fn test_extract_against_fasta_fixtures() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_fasta = temp_dir.path().join("out.fasta");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the extract command
        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/simple.fasta"),
            in_fastq_2: None,
            kmer_seq: Some(vec!["ACG".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_fasta.clone()),
            q_size: None,
            aho_corasick: false,
            reverse_complement: true,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 1,
        };

        extract_records(args)?;

        // Compare outputs with fixtures
        compare_fasta_output(&out_fasta, "tests/fixtures/extract/simple.extracted.fasta")?;
        compare_log_output(&out_log, "tests/fixtures/extract/simple.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/simple.json")?;

        Ok(())
    }

    // Test inverted matching mode with simple nucleotide FASTA file
    // Corresponds to: cargo run -- extract -i tests/fixtures/input/simple.fasta -r -s ACG -v -o tests/fixtures/extract/simple-inv.extracted.fasta -l tests/fixtures/extract/simple-inv.log -j tests/fixtures/extract/simple-inv.json
    #[test]
    fn test_extract_against_fasta_fixtures_inverted() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_fasta = temp_dir.path().join("out.fasta");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the extract command with inverted matching
        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/simple.fasta"),
            in_fastq_2: None,
            kmer_seq: Some(vec!["ACG".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_fasta.clone()),
            q_size: None,
            aho_corasick: false,
            reverse_complement: true,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: true,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 1,
        };

        extract_records(args)?;

        // Compare outputs with fixtures
        compare_fasta_output(
            &out_fasta,
            "tests/fixtures/extract/simple-inv.extracted.fasta",
        )?;
        compare_log_output(&out_log, "tests/fixtures/extract/simple-inv.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/simple-inv.json")?;

        Ok(())
    }

    // Compare with fixed-width amino acid FASTA file, containing a match at a line break
    // Corresponds to: cargo run -- extract -i tests/fixtures/input/fixed-width.faa -s DKAT -o tests/fixtures/extract/fixed-width.extracted.faa -l tests/fixtures/extract/fixed-width.log -j tests/fixtures/extract/fixed-width.json
    #[test]
    fn test_extract_against_fasta_fixtures_fixed_width_aa() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_fasta = temp_dir.path().join("out.faa");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the extract command
        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/fixed-width.faa"),
            in_fastq_2: None,
            kmer_seq: Some(vec!["DKAT".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_fasta.clone()),
            q_size: None,
            aho_corasick: false,
            reverse_complement: false,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 1,
        };

        extract_records(args)?;

        // Compare outputs with fixtures
        compare_fasta_output(
            &out_fasta,
            "tests/fixtures/extract/fixed-width.extracted.faa",
        )?;
        compare_log_output(&out_log, "tests/fixtures/extract/fixed-width.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/fixed-width.json")?;

        Ok(())
    }

    // Compare with paired-end FASTQ files
    // Corresponds to: cargo run -- extract -i tests/fixtures/input/paired-1.fastq -2 tests/fixtures/input/paired-2.fastq -s CTT -o tests/fixtures/extract/paired.extracted.fastq -l tests/fixtures/extract/paired.log -j tests/fixtures/extract/paired.json
    #[test]
    fn test_extract_against_fastq_fixtures_paired() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_base = temp_dir.path().join("out");
        let out_fastq_1 = temp_dir.path().join("out_1.fastq");
        let out_fastq_2 = temp_dir.path().join("out_2.fastq");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the extract command
        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/paired-1.fastq"),
            in_fastq_2: Some(PathBuf::from("tests/fixtures/input/paired-2.fastq")),
            kmer_seq: Some(vec!["CTT".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_base),
            q_size: None,
            aho_corasick: false,
            reverse_complement: false,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 1,
        };

        extract_records(args)?;

        // Compare outputs with fixtures
        compare_fasta_output(
            &out_fastq_1,
            "tests/fixtures/extract/paired_1.extracted.fastq",
        )?;
        compare_fasta_output(
            &out_fastq_2,
            "tests/fixtures/extract/paired_2.extracted.fastq",
        )?;
        compare_log_output(&out_log, "tests/fixtures/extract/paired.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/paired.json")?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_single_end_matches_serial_fixtures() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let out_fasta = temp_dir.path().join("out.fasta");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/simple.fasta"),
            in_fastq_2: None,
            kmer_seq: Some(vec!["ACG".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_fasta.clone()),
            q_size: None,
            aho_corasick: false,
            reverse_complement: true,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 4,
        };

        extract_records(args)?;

        compare_fasta_output(&out_fasta, "tests/fixtures/extract/simple.extracted.fasta")?;
        compare_log_output(&out_log, "tests/fixtures/extract/simple.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/simple.json")?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_paired_end_matches_serial_fixtures() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let out_base = temp_dir.path().join("out");
        let out_fastq_1 = temp_dir.path().join("out_1.fastq");
        let out_fastq_2 = temp_dir.path().join("out_2.fastq");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/paired-1.fastq"),
            in_fastq_2: Some(PathBuf::from("tests/fixtures/input/paired-2.fastq")),
            kmer_seq: Some(vec!["CTT".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_base),
            q_size: None,
            aho_corasick: false,
            reverse_complement: false,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 4,
        };

        extract_records(args)?;

        compare_fasta_output(
            &out_fastq_1,
            "tests/fixtures/extract/paired_1.extracted.fastq",
        )?;
        compare_fasta_output(
            &out_fastq_2,
            "tests/fixtures/extract/paired_2.extracted.fastq",
        )?;
        compare_log_output(&out_log, "tests/fixtures/extract/paired.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/paired.json")?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_threads_2_and_4_match_serial_single_end() -> Result<()> {
        let serial = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 1))?;
        let parallel_2 =
            run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 2))?;
        let parallel_4 =
            run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;

        assert_single_runs_match(&parallel_2, &serial)?;
        assert_single_runs_match(&parallel_4, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_threads_4_is_deterministic() -> Result<()> {
        let first = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;
        let second = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;
        let third = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;

        assert_single_runs_match(&second, &first)?;
        assert_single_runs_match(&third, &first)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_inverted_matches_serial() -> Result<()> {
        let mut serial_options = SingleExtractOptions::simple(vec!["ACG".to_string()], 1);
        serial_options.invert_match = true;
        let mut parallel_options = serial_options.clone();
        parallel_options.threads = 4;

        let serial = run_single_extract(serial_options)?;
        let parallel = run_single_extract(parallel_options)?;

        assert_single_runs_match(&parallel, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_aho_corasick_matches_serial() -> Result<()> {
        let mut serial_options = SingleExtractOptions::simple(vec!["ACG".to_string()], 1);
        serial_options.aho_corasick = true;
        let mut parallel_options = serial_options.clone();
        parallel_options.threads = 4;

        let serial = run_single_extract(serial_options)?;
        let parallel = run_single_extract(parallel_options)?;

        assert_single_runs_match(&parallel, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_bndmq_matches_serial() -> Result<()> {
        let serial = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 1))?;
        let parallel =
            run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;

        assert_single_runs_match(&parallel, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_case_insensitive_matches_serial() -> Result<()> {
        let mut serial_options = SingleExtractOptions::simple(vec!["acg".to_string()], 1);
        serial_options.case_insensitive = true;
        let mut parallel_options = serial_options.clone();
        parallel_options.threads = 4;

        let serial = run_single_extract(serial_options)?;
        let parallel = run_single_extract(parallel_options)?;

        assert_single_runs_match(&parallel, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_fixed_width_fasta_matches_serial() -> Result<()> {
        let mut serial_options = SingleExtractOptions::simple(vec!["DKAT".to_string()], 1);
        serial_options.input = PathBuf::from("tests/fixtures/input/fixed-width.faa");
        serial_options.reverse_complement = false;
        let mut parallel_options = serial_options.clone();
        parallel_options.threads = 4;

        let serial = run_single_extract(serial_options)?;
        let parallel = run_single_extract(parallel_options)?;

        assert_single_runs_match(&parallel, &serial)?;

        Ok(())
    }

    #[test]
    fn test_extract_parallel_no_matches() -> Result<()> {
        let mut options = SingleExtractOptions::simple(vec!["GGG".to_string()], 4);
        options.reverse_complement = false;
        let run = run_single_extract(options)?;
        let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&run.out_json)?)?;

        assert_eq!(fs::read_to_string(&run.out_fastx)?, "");
        assert_eq!(json["summary_statistics"]["number_of_matches"], 0);
        assert_eq!(
            json["paired_end_reads_statistics"]["number_of_extracted_records"],
            0
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_all_records_match() -> Result<()> {
        let mut options = SingleExtractOptions::simple(vec!["T".to_string()], 4);
        options.reverse_complement = false;
        let run = run_single_extract(options)?;
        let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&run.out_json)?)?;

        assert_eq!(json["summary_statistics"]["number_of_records_searched"], 3);
        assert_eq!(
            json["paired_end_reads_statistics"]["number_of_extracted_records"],
            3
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_dense_matches_with_verbose_logging() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let input = temp_dir.path().join("dense.fasta");
        fs::write(&input, b">dense\nAAAAA\n>sparse\nTTTTA\n")?;

        let mut options = SingleExtractOptions::simple(vec!["A".to_string()], 4);
        options.input = input;
        options.reverse_complement = false;
        let run = run_single_extract(options)?;
        let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&run.out_json)?)?;

        assert_eq!(json["summary_statistics"]["number_of_records_searched"], 2);
        assert_eq!(json["summary_statistics"]["number_of_matches"], 6);
        assert_eq!(
            json["paired_end_reads_statistics"]["number_of_extracted_records"],
            2
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_sparse_matches_with_verbose_logging() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let input = temp_dir.path().join("sparse.fasta");
        fs::write(&input, b">hit\nACGT\n>miss1\nTTTT\n>miss2\nCCCC\n")?;

        let mut options = SingleExtractOptions::simple(vec!["ACG".to_string()], 4);
        options.input = input;
        options.reverse_complement = false;
        let run = run_single_extract(options)?;
        let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&run.out_json)?)?;

        assert_eq!(json["summary_statistics"]["number_of_records_searched"], 3);
        assert_eq!(json["summary_statistics"]["number_of_matches"], 1);
        assert_eq!(
            json["paired_end_reads_statistics"]["number_of_extracted_records"],
            1
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_suppress_output_logs_only() -> Result<()> {
        let mut options = SingleExtractOptions::simple(vec!["ACG".to_string()], 4);
        options.suppress_output = true;
        let run = run_single_extract(options)?;
        let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(&run.out_json)?)?;

        assert!(!run.out_fastx.exists());
        assert!(run.out_log.exists());
        assert_eq!(
            json["paired_end_reads_statistics"]["number_of_extracted_records"],
            2
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_compressed_fasta_inputs_match_uncompressed() -> Result<()> {
        let mut uncompressed_options = SingleExtractOptions::simple(vec!["ATGG".to_string()], 4);
        uncompressed_options.input = PathBuf::from("tests/data/sample.fasta");
        uncompressed_options.reverse_complement = false;
        let uncompressed = run_single_extract(uncompressed_options.clone())?;

        for compressed_path in [
            "tests/data/sample.fasta.gz",
            "tests/data/sample.fasta.bz2",
            "tests/data/sample.fasta.xz",
        ] {
            let mut compressed_options = uncompressed_options.clone();
            compressed_options.input = PathBuf::from(compressed_path);
            let compressed = run_single_extract(compressed_options)?;
            compare_text_files(&compressed.out_fastx, &uncompressed.out_fastx)?;

            let expected_json: serde_json::Value =
                serde_json::from_str(&fs::read_to_string(&uncompressed.out_json)?)?;
            let actual_json: serde_json::Value =
                serde_json::from_str(&fs::read_to_string(&compressed.out_json)?)?;
            assert_eq!(
                expected_json["summary_statistics"],
                actual_json["summary_statistics"]
            );
            assert_eq!(
                expected_json["paired_end_reads_statistics"],
                actual_json["paired_end_reads_statistics"]
            );
            assert_eq!(
                expected_json["pattern_hit_counts"],
                actual_json["pattern_hit_counts"]
            );
        }

        Ok(())
    }

    #[test]
    fn test_extract_thread_resolution_auto_detects_available_threads() {
        assert_eq!(
            resolve_extract_threads_with_available(0, 8),
            ThreadResolution {
                effective_threads: 8,
                auto_detected: true,
                clamped: false,
            }
        );
    }

    #[test]
    fn test_extract_thread_resolution_clamps_to_available_threads() {
        assert_eq!(
            resolve_extract_threads_with_available(16, 8),
            ThreadResolution {
                effective_threads: 8,
                auto_detected: false,
                clamped: true,
            }
        );
    }

    #[test]
    fn test_extract_auto_threads_matches_serial_fixtures() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let out_fasta = temp_dir.path().join("out.fasta");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        let args = CmdExtract {
            in_fastx: PathBuf::from("tests/fixtures/input/simple.fasta"),
            in_fastq_2: None,
            kmer_seq: Some(vec!["ACG".to_string()]),
            kmer_file: None,
            out_fastx: Some(out_fasta.clone()),
            q_size: None,
            aho_corasick: false,
            reverse_complement: true,
            canonical: false,
            out_log: Some(out_log.clone()),
            suppress_output: false,
            json_log: Some(out_json.clone()),
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 0,
        };

        extract_records(args)?;

        compare_fasta_output(&out_fasta, "tests/fixtures/extract/simple.extracted.fasta")?;
        compare_log_output(&out_log, "tests/fixtures/extract/simple.log")?;
        compare_json_output(&out_json, "tests/fixtures/extract/simple.json")?;

        Ok(())
    }

    #[test]
    fn test_extract_paired_fastq_length_mismatch_errors() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let read_1 = temp_dir.path().join("read_1.fastq");
        let read_2 = temp_dir.path().join("read_2.fastq");
        fs::write(
            &read_1,
            b"@seq1/1\nACTTACGT\n+\nIIIIIIII\n@seq2/1\nTTTTTTTT\n+\nIIIIIIII\n",
        )?;
        fs::write(&read_2, b"@seq1/2\nGCTATAAT\n+\nIIIIIIII\n")?;

        let args = CmdExtract {
            in_fastx: read_1,
            in_fastq_2: Some(read_2),
            kmer_seq: Some(vec!["CTT".to_string()]),
            kmer_file: None,
            out_fastx: None,
            q_size: None,
            aho_corasick: false,
            reverse_complement: false,
            canonical: false,
            out_log: None,
            suppress_output: false,
            json_log: None,
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 1,
        };

        let error = extract_records(args).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("The two input files have a different number of records")
        );

        Ok(())
    }

    #[test]
    fn test_extract_parallel_paired_fastq_length_mismatch_errors() -> Result<()> {
        let temp_dir = tempfile::tempdir()?;
        let read_1 = temp_dir.path().join("read_1.fastq");
        let read_2 = temp_dir.path().join("read_2.fastq");
        fs::write(
            &read_1,
            b"@seq1/1\nACTTACGT\n+\nIIIIIIII\n@seq2/1\nTTTTTTTT\n+\nIIIIIIII\n",
        )?;
        fs::write(&read_2, b"@seq1/2\nGCTATAAT\n+\nIIIIIIII\n")?;

        let args = CmdExtract {
            in_fastx: read_1,
            in_fastq_2: Some(read_2),
            kmer_seq: Some(vec!["CTT".to_string()]),
            kmer_file: None,
            out_fastx: None,
            q_size: None,
            aho_corasick: false,
            reverse_complement: false,
            canonical: false,
            out_log: None,
            suppress_output: false,
            json_log: None,
            invert_match: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
            threads: 4,
        };

        let error = extract_records(args).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("The two input files have a different number of records")
        );

        Ok(())
    }
}
