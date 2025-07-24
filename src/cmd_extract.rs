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

use aho_corasick::AhoCorasick;
use anyhow::{Context, Result};
use clap::{ArgAction, ArgGroup, Args, crate_name, crate_version};
use jiff::{Unit, Zoned};
use serde_json;

use std::collections::HashMap;
use std::{fs, env};
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::string::String;

use crate::helpers::{
    add_suffix_to_file_prefix, check_log_flag_conflict, identify_uncompressed_type,
    parse_pattern_list, recommend_aho_corasick, error_if_directory,
};
use crate::logger::{BufferedLogger, JsonLogger};
use crate::pattern_matching::{BNDMq, tune_q_value};

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
#[derive(Clone)]
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
        },
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
        logger.write_header("#SeqKatcher extract log\n");
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

    // Initialize algorithm instances for each pattern. Only construct the Aho-
    // Corasick automaton when requested.
    let (ac, bndmq_collection): (Option<AhoCorasick>, Vec<(String, BNDMq)>) = if args.aho_corasick {
        let ac = AhoCorasick::builder()
            // Use DFA for better search performance at higher memory cost
            .kind(Some(aho_corasick::AhoCorasickKind::DFA))
            .ascii_case_insensitive(args.case_insensitive)
            .build(pattern_list.clone())
            .unwrap();
        (Some(ac), Vec::new())
    } else {
        let mut bndmq_collection = Vec::with_capacity(pattern_list.len());
        for pattern in &pattern_list {
            // Tune q value for BNDMq based on the pattern length if not provided
            let q = args
                .q_size
                .unwrap_or_else(|| tune_q_value(pattern).unwrap());
            bndmq_collection.push((pattern.clone(), BNDMq::new(pattern.as_bytes(), q)?));
        }
        (None, bndmq_collection)
    };

    // Uses a gzip decoder or regular file reader to read FASTQ/A records,
    // depending on the file extension
    let mut reader = needletail::parse_fastx_file(&args.in_fastx)
        .with_context(|| format!("Invalid FASTQ/A input path or file: {:?}", &args.in_fastx))?;

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
        // Either write to file or stdout if no output path is provided;
        // the file format is determined by the input file
        let mut writer = match &args.out_fastx {
            Some(pathbuf) => {
                let pathbuf =
                    pathbuf.with_extension(identify_uncompressed_type(&args.in_fastx).unwrap());
                let path = Path::new(&pathbuf);
                let file = fs::File::create(path).with_context(|| {
                    format!(
                        "Error writing to output file; no such directory: {path:?}",
                        
                    )
                })?;
                Box::new(BufWriter::new(file)) as Box<dyn io::Write>
            }
            None => {
                let stdout = io::stdout();
                Box::new(BufWriter::new(stdout)) as Box<dyn io::Write>
            }
        };

        // Iterate over FASTA/Q records and check for k-mer presence
        while let Some(r) = reader.next() {
            let record = r.with_context(|| "Error during FASTQ/A record parsing.")?;
            let mut found_occ = false;

            if logging_active {
                nb_records_tot += 1;
                nb_bases += record.num_bases();
            }

            // Get occurrences of k-mers in the sequence using Aho-Corasick
            if let Some(ac) = ac.as_ref() {
                for mat in ac.find_overlapping_iter(&record.seq()) {
                    if !logging_active {
                        found_occ = true;
                        break;
                    } else {
                        if logging_active {
                            logger.log_fields(
                                in_fastx_filename,
                                record.id(),
                                &pattern_list[mat.pattern().as_usize()],
                                mat.start(),
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastx_filename,
                                    record.id(),
                                    &pattern_list[mat.pattern().as_usize()],
                                    mat.start(),
                                );
                            }
                        }
                        pattern_hit_counts[mat.pattern().as_usize()] += 1;
                        nb_hits_tot.0 += 1;
                        found_occ = true;
                    }
                }
                if found_occ {
                    nb_records_hit.0 += 1;
                }
            // Or use BNDMq
            } else {
                // If logging active, search for matching positions and print them
                if logging_active {
                    for (idx, (pattern, bndmq)) in bndmq_collection.iter().enumerate() {
                        let mut found_any = false;
                        for o in bndmq.find_iter(&record.seq()) {
                            found_any = true;
                            logger.log_fields(
                                in_fastx_filename,
                                record.id(),
                                pattern,
                                o,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(in_fastx_filename, record.id(), pattern, o);
                            }
                            nb_hits_tot.0 += 1;
                        }
                        if found_any {
                            found_occ = true;
                            pattern_hit_counts[idx] += 1;
                        }
                    }
                    if found_occ {
                        nb_records_hit.0 += 1;
                    }
                // If logging disabled, only search for a match and break if found
                } else {
                    for (_, bndmq) in &bndmq_collection {
                        if bndmq.find_match(&record.seq()) {
                            found_occ = true;
                            break;
                        }
                    }
                }
            }

            // Write record to file or stdout if any k-mer has been found
            if found_occ != args.invert_match {
                nb_records_extracted += 1;
                if !args.suppress_output {
                    record.write(&mut writer, None).unwrap();
                }
            }
        }
    ////
    //// ---------------------- Handling Paired-End Reads ---------------------
    ////
    // If a second file is provided, process paired-end reads
    } else {
        let mut reader_2 = needletail::parse_fastx_file(args.in_fastq_2.clone().unwrap())
            .with_context(|| {
                format!(
                    "Invalid second FASTQ input path or file: {:?}",
                    &args.in_fastq_2
                )
            })?;

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
                    format!(
                        "Error writing to paired-end file; no such directory: {path:?}",
                        
                    )
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
                    format!(
                        "Error writing second paired-end file; no such directory: {path:?}",
                        
                    )
                })?;
                Box::new(BufWriter::new(file)) as Box<dyn io::Write>
            }
            None => {
                let stdout = io::stdout();
                Box::new(BufWriter::new(stdout)) as Box<dyn io::Write>
            }
        };

        // Iterate over FASTQ records and check for k-mer presence
        while let Some(r) = reader.next() {
            let record_1 = r.with_context(|| "Error during FASTQ record parsing of first file.")?;
            let record_2 = reader_2
                .next()
                .with_context(|| "Error during FASTQ record parsing of second file. Do the two input files contain the same number of records?")?
                .unwrap();
            let mut found_occ = false;

            if logging_active {
                nb_records_tot += 2;
                nb_bases += record_1.num_bases();
                nb_bases += record_2.num_bases();
            }

            // Get occurrences of patterns in the sequence using Aho-Corasick
            if let Some(ac) = ac.as_ref() {
                let mut record_hit: (usize, usize) = (0, 0);
                for mat in ac.find_overlapping_iter(&record_1.seq()) {
                    if !logging_active {
                        found_occ = true;
                        break;
                    } else {
                        if logging_active {
                            logger.log_fields(
                                in_fastx_filename,
                                record_1.id(),
                                &pattern_list[mat.pattern().as_usize()],
                                mat.start(),
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastx_filename,
                                    record_1.id(),
                                    &pattern_list[mat.pattern().as_usize()],
                                    mat.start(),
                                );
                            }
                        }
                        pattern_hit_counts[mat.pattern().as_usize()] += 1;
                        record_hit.0 = 1;
                        nb_hits_tot.0 += 1;
                        found_occ = true;
                    }
                }
                for mat in ac.find_overlapping_iter(&record_2.seq()) {
                    if !logging_active {
                        found_occ = true;
                        break;
                    } else {
                        if logging_active {
                            logger.log_fields(
                                in_fastq_2_filename,
                                record_2.id(),
                                &pattern_list[mat.pattern().as_usize()],
                                mat.start(),
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(
                                    in_fastq_2_filename,
                                    record_2.id(),
                                    &pattern_list[mat.pattern().as_usize()],
                                    mat.start(),
                                );
                            }
                        }
                        pattern_hit_counts[mat.pattern().as_usize()] += 1;
                        record_hit.1 = 1;
                        nb_hits_tot.1 += 1;
                        found_occ = true;
                    }
                }
                if logging_active {
                    nb_records_hit.0 += record_hit.0;
                    nb_records_hit.1 += record_hit.1;
                }
            // Or use BNDMq
            } else {
                // If logging active, search for matching positions and print them
                if logging_active {
                    let mut record_hit: (usize, usize) = (0, 0);
                    for (idx, (pattern, bndmq)) in bndmq_collection.iter().enumerate() {
                        let mut found_any1 = false;
                        let mut found_any2 = false;

                        for o in bndmq.find_iter(&record_1.seq()) {
                            found_any1 = true;
                            logger.log_fields(
                                in_fastx_filename,
                                record_1.id(),
                                pattern,
                                o,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(in_fastx_filename, record_1.id(), pattern, o);
                            }
                            nb_hits_tot.0 += 1;
                        }

                        for o in bndmq.find_iter(&record_2.seq()) {
                            found_any2 = true;
                            logger.log_fields(
                                in_fastq_2_filename,
                                record_2.id(),
                                pattern,
                                o,
                            );
                            if let Some(jl) = &mut json_logger {
                                jl.log_fields(in_fastq_2_filename, record_2.id(), pattern, o);
                            }
                            nb_hits_tot.1 += 1;
                        }

                        if found_any1 {
                            found_occ = true;
                            record_hit.0 = 1;
                            pattern_hit_counts[idx] += 1;
                        }
                        if found_any2 {
                            found_occ = true;
                            record_hit.1 = 1;
                            pattern_hit_counts[idx] += 1;
                        }
                    }
                    nb_records_hit.0 += record_hit.0;
                    nb_records_hit.1 += record_hit.1;
                // If logging disabled, only search for a match and break if found
                } else {
                    for (_, bndmq) in &bndmq_collection {
                        if bndmq.find_match(&record_1.seq()) || bndmq.find_match(&record_2.seq()) {
                            found_occ = true;
                            break;
                        }
                    }
                }
            }

            // Write records to file or stdout if any patterns have been matched
            if found_occ != args.invert_match {
                nb_records_extracted += 2;
                if !args.suppress_output {
                    record_1.write(&mut writer, None).unwrap();
                    record_2.write(&mut writer2, None).unwrap();
                }
            }
        }
        if reader_2.next().is_some() {
            anyhow::bail!(
                "The two input files have a different number of records. Please provide valid paired-end read files."
            );
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
        let pattern_hit_counts_map: HashMap<String, u32> =
            pattern_list.iter().cloned().zip(pattern_hit_counts.iter().copied()).collect();
        let meta_information = serde_json::json!({
            "program": crate_name!(),
            "version": crate_version!(),
            "timestamp": Zoned::now().round(Unit::Second).unwrap(),
            "subcommand": "extract",
            "command_line": env::args().collect::<Vec<String>>(),
            "search_algorithm": if args.aho_corasick { "Aho-Corasick" } else { "BNDMq" },
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

    /// Compare FASTA output with fixture
    fn compare_fasta_output(actual_path: &Path, expected_path: &str) -> Result<()> {
        let expected = fs::read_to_string(expected_path)?;
        let actual = fs::read_to_string(actual_path)?;
        assert_eq!(expected, actual, "FASTA output does not match fixture");
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
}
