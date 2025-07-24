//! # Tag subcommand
//!
//! The `tag` subcommand reads a BAM or SAM file, searches for subsequences
//! and tags the records with the presence of k-mers. The output is written
//! to a new BAM or SAM file, adding a tag to records. File type is determined
//! automatically based on the input file extension. Additionally, adds a tag to the BAM/SAM
//! header with the program information.

use aho_corasick::AhoCorasick;
use anyhow::{Context, Result};
use bam::record::tags;
use bam::{RecordReader, RecordWriter};
use clap::{ArgAction, ArgGroup, Args, crate_name, crate_version};
use jiff::{Unit, Zoned};
use serde_json;

use std::collections::HashMap;
use std::path::PathBuf;
use std::str::from_utf8;
use std::{env, fs, io};

use crate::helpers::{
    check_log_flag_conflict, error_if_directory, parse_pattern_list, recommend_aho_corasick,
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
    ArgGroup::new("logging")
        .required(false)
        .multiple(true)
        .args(&["out_log", "json_log"]),
),
group(
    ArgGroup::new("matching")
        .required(false)
        .multiple(false)
        .args(&["filter_matching", "invert_match"]),
),
group(
    ArgGroup::new("algorithm")
        .required(false)
        .multiple(false)
        .args(&["q_size", "aho_corasick"]),
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
pub struct CmdTag {
    /// Input path for SAM/BAM file.
    #[clap(short = 'i', long)]
    in_file: PathBuf,

    /// Output path for SAM/BAM file with annotations; file type is inferred from the file extension.
    #[clap(short = 'o', long, required = false)]
    out_file: Option<PathBuf>,

    /// Query sequences (accepts multiple sequences after the flag, separated by a space); if not provided, input path for file containing list of k-mers is required.
    #[clap(short = 's', long, num_args = 1..)]
    kmer_seq: Option<Vec<String>>,

    /// Input path for file containing list of k-mers, one per line (FASTA or plain text file; comment lines starting with '#'are ignored).
    #[clap(short = 'f', long)]
    kmer_file: Option<PathBuf>,

    /// Also search for reverse complements of k-mers.
    #[clap(short = 'r', long, action(ArgAction::SetTrue), default_value("false"))]
    reverse_complement: bool,

    /// Search only for the canonical forms of k-mers.
    #[clap(short = 'c', long, action(ArgAction::SetTrue), default_value("false"))]
    canonical: bool,

    /// Tag to add to the SAM/BAM file with the presence of k-mers.
    #[clap(short = 't', long, default_value("km"))]
    tag: String,

    /// Print detailed match information to stdout if only the flag is passed, or to a file if a path is provided.
    #[clap(short = 'l', long, default_value(None), default_missing_value("STDOUT"), num_args = 0..=1, )]
    out_log: Option<PathBuf>,

    /// Write JSON log to stdout if only the flag is passed, or to a file if a path is provided.
    #[clap(short = 'j', long, default_value(None), default_missing_value("STDOUT"), num_args = 0..=1, )]
    json_log: Option<PathBuf>,

    /// Number of parallel threads to use for processing BAM files.
    #[clap(short = 'p', long, default_value("1"))]
    threads: u16,

    /// Suppress output of found records (no records are written to a file or stdout); use if only matching statistics are of interest.
    #[clap(
        short = 'S',
        long,
        action(ArgAction::SetTrue),
        default_value("false"),
        conflicts_with("out_file"),
        requires("logging")
    )]
    suppress_output: bool,

    /// Filter records to keep only those with matching k-mers.
    #[clap(short = 'm', long, action(ArgAction::SetTrue), default_value("false"))]
    filter_matching: bool,

    /// Invert the sense of matching, filtering out records that match instead of keeping them.
    #[clap(short = 'v', long, action(ArgAction::SetTrue), default_value("false"))]
    invert_match: bool,

    /// Use case-insensitive matching.
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

/// Core function of the `tag` subcommand that reads a SAM/BAM file, searches
/// for subsequences and tags the records with the presence of k-mers.
/// The output is written to a new SAM/BAM file, adding a tag to records.
pub fn tag_records(args: CmdTag) -> Result<()> {
    // Use helper for log flag conflict (not possible yet with `clap`)
    check_log_flag_conflict(
        &args.out_log,
        &args.json_log,
        &args.out_file,
        args.suppress_output,
    )
    .map_err(|e| anyhow::anyhow!(e))?;

    let mut args = args;

    // Check if input file path points to a directory
    error_if_directory(&args.in_file, "Record file path")?;

    // Get input filename for logging
    let in_records_filename = args.in_file.file_name().unwrap().to_str().unwrap();

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
    let log_file =
        match &args.out_log {
            Some(path) => {
                if path
                    .to_str()
                    .ok_or_else(|| anyhow::anyhow!("Invalid log file path."))?
                    == "STDOUT"
                {
                    Some(Box::new(io::stdout()) as Box<dyn io::Write>)
                } else {
                    Some(Box::new(fs::File::create(path).with_context(|| {
                        format!("Problem creating log file: {}", path.display())
                    })?) as Box<dyn io::Write>)
                }
            }
            None => None,
        };

    let out_file = args.out_file.clone();

    // Activate logging if a log or JSON log file is provided
    let logging_active = log_file.is_some() || args.json_log.is_some();

    // Check if number of threads is at least 1
    if args.threads < 1 {
        anyhow::bail!("Number of threads must be at least 1.");
    }
    // Check if the tag is a valid tag name
    let tag_validated: [u8; 2] = if args.tag.len() != 2 {
        anyhow::bail!("Tag must be exactly two characters long.");
    } else {
        args.tag
            .as_bytes()
            .try_into()
            .map_err(|_| anyhow::anyhow!("Invalid tag format."))?
    };

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

    fn infer_record_writer(
        threads: u16,
        out_file: &Option<PathBuf>,
        extension: &str,
        header: bam::Header,
    ) -> Result<Box<dyn RecordWriter>> {
        match extension {
            "bam" => {
                let out_file = out_file
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("Output file not provided for BAM writing."))?;
                let path = out_file.with_extension(extension);
                Ok(Box::new(
                    bam::bam_writer::BamWriterBuilder::new()
                        .additional_threads(threads - 1)
                        .from_path(&path, header)
                        .with_context(|| format!("Error writing BAM file: {}", path.display()))?,
                ))
            }
            "sam" => {
                let out_file = out_file
                    .as_ref()
                    .ok_or_else(|| anyhow::anyhow!("Output file not provided for SAM writing."))?;
                let path = out_file.with_extension(extension);
                Ok(Box::new(
                    bam::sam::SamWriterBuilder::new()
                        .from_path(&path, header)
                        .with_context(|| format!("Error writing SAM file: {}", path.display()))?,
                ))
            }
            "STDOUT" => Ok(Box::new(
                bam::sam::SamWriterBuilder::new()
                    .from_stream(io::stdout(), header)
                    .with_context(|| "Error writing SAM file to stdout.")?,
            )),
            _ => anyhow::bail!("Invalid output file type. Must be BAM or SAM file."),
        }
    }

    // Check if file is a BAM or SAM file and open it for reading
    let in_file_extension = args
        .in_file
        .extension()
        .with_context(|| format!("Could not detect the file extension: {:?}", args.in_file))?
        .to_str()
        .unwrap();

    // If out_file is provided, set depending on extension; otherwise, set to STDOUT
    let out_file_extension = match &args.out_file {
        Some(path) => match path.extension() {
            Some(extension) => extension.to_str().unwrap(),
            None => in_file_extension,
        },
        None => "STDOUT",
    };

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
        logger.write_header("#SeqKatcher tag log\n");
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
            "#Tag used for labeling records: {}\n",
            from_utf8(&tag_validated).unwrap()
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

    // Initialize counters for logging information
    let mut nb_records_tot = 0;
    let mut nb_bases: usize = 0;
    let mut nb_hits_tot = 0;
    let mut nb_records_hit = 0;
    let mut pattern_hit_counts = vec![0u32; pattern_list.len()];

    /// Process a single record, updating statistics and writing to output
    fn process_record(
        record: &mut bam::Record,
        ac: Option<&AhoCorasick>,
        bndmq_collection: &[(String, BNDMq)],
        pattern_list: &[String],
        tag_validated: [u8; 2],
        logging_active: bool,
        logger: &mut BufferedLogger,
        in_records_filename: &str,
        pattern_hit_counts: &mut [u32],
        nb_hits_tot: &mut usize,
        nb_records_hit: &mut usize,
        nb_records_tot: &mut usize,
        nb_bases: &mut usize,
        filter_matching: bool,
        suppress_output: bool,
        invert_match: bool,
        writer: &mut Box<dyn RecordWriter>,
        json_logger: &mut Option<JsonLogger>,
    ) -> Result<()> {
        let mut kmers_found = Vec::new();

        let use_ac = ac.is_some();

        // Get occurrences of patterns in the sequence using Aho-Corasick or BNDMq
        if use_ac {
            for mat in ac
                .unwrap()
                .find_overlapping_iter(&record.sequence().to_vec())
            {
                if let Some(pattern) = pattern_list.get(mat.pattern().as_usize()) {
                    kmers_found.push(pattern.clone());
                    // Log match information
                    if logging_active {
                        *nb_hits_tot += 1;

                        if let Some(count) = pattern_hit_counts.get_mut(mat.pattern().as_usize()) {
                            *count += 1;
                        }
                        logger.log_fields(in_records_filename, record.name(), pattern, mat.start());
                        if let Some(jl) = json_logger.as_mut() {
                            jl.log_fields(in_records_filename, record.name(), pattern, mat.start());
                        }
                    }
                } else {
                    anyhow::bail!("Error retrieving matching pattern by index.")
                }
            }
        } else {
            // If logging active, search for matching positions and print them
            if logging_active {
                for (idx, (pattern, bndmq)) in bndmq_collection.iter().enumerate() {
                    let mut found_any = false;
                    for o in bndmq.find_iter(&record.sequence().to_vec()) {
                        found_any = true;
                        logger.log_fields(in_records_filename, record.name(), pattern, o);
                        if let Some(jl) = json_logger.as_mut() {
                            jl.log_fields(in_records_filename, record.name(), pattern, o);
                        }
                        *nb_hits_tot += 1;
                    }
                    if found_any {
                        kmers_found.push(pattern.clone());
                        if let Some(count) = pattern_hit_counts.get_mut(idx) {
                            *count += 1;
                        }
                    }
                }
            // If logging disabled, only search for a match and break if found
            } else {
                for (pattern, bndmq) in bndmq_collection {
                    if bndmq.find_match(&record.sequence().to_vec()) {
                        kmers_found.push(pattern.clone());
                    }
                }
            }
        }

        if logging_active {
            *nb_records_tot += 1;
            *nb_bases += record.query_len() as usize;
            if !kmers_found.is_empty() {
                *nb_records_hit += 1;
            }
        }

        // Skip record based on matching criteria:
        // - With filter_matching (-m): keep only records that match
        // - With invert_match (-v): keep only records that don't match
        // - Without either: keep all records
        let should_keep = if filter_matching {
            !kmers_found.is_empty() // Keep only matching records
        } else if invert_match {
            kmers_found.is_empty() // Keep only non-matching records
        } else {
            true // Keep all records
        };

        if !should_keep {
            return Ok(());
        }

        // Tag record with presence of k-mers
        match record.tags().get(&tag_validated) {
            // Do nothing if tag is empty
            Some(tags::TagValue::String([], _)) => (),
            // Otherwise, append the new k-mers to the newly found k-mers
            Some(tags::TagValue::String(val, _)) => {
                let s =
                    from_utf8(val).with_context(|| "Error reading existing tag value as UTF-8")?;
                kmers_found.extend(s.split(',').map(String::from));
            }
            None => (),
            _ => anyhow::bail!("Invalid tag value format. Expected string value."),
        };

        // Sort and deduplicate k-mers
        kmers_found.sort_unstable();
        kmers_found.dedup();

        // Update record with new k-mers
        record
            .tags_mut()
            .push_string(&tag_validated, kmers_found.join(",").as_bytes());

        // Write record to output file if not suppressed
        if !suppress_output {
            writer
                .write(record)
                .with_context(|| "Error writing record to output file")?;
        }

        Ok(())
    }

    // Process records based on file type
    match in_file_extension {
        "bam" => {
            // Open BAM file for reading with x additional threads for decompression
            let mut reader = bam::BamReader::from_path(&args.in_file, &args.threads - 1)
                .with_context(|| format!("Error reading BAM file: {:?}", &args.in_file))?;
            // Get header from BAM file and add program information
            let command_line = env::args().collect::<Vec<String>>().join(" ");
            let mut program_header_line = format!("@PG\tID:{0}\tPN:{0}\tCL:", crate_name!());
            program_header_line.push_str(&command_line);
            program_header_line.push_str(format!("\tVN:{}", crate_version!()).as_str());
            let mut header = reader.header().clone();
            header.push_line(&program_header_line).unwrap();
            // Use empty header if suppress_output is set
            if args.suppress_output {
                header = bam::Header::new();
            }
            // Open file for writing with inferred writer
            let mut writer = match out_file_extension {
                "bam" | "sam" | "STDOUT" => {
                    infer_record_writer(args.threads, &out_file, out_file_extension, header)
                }
                _ => anyhow::bail!("Output file must be a BAM or SAM file."),
            }
            .with_context(|| "Could not create writer.")?;

            // Iterate over BAM records and process them
            let mut record = bam::Record::new();
            loop {
                match reader.read_into(&mut record) {
                    Ok(true) => {
                        process_record(
                            &mut record,
                            ac.as_ref(),
                            &bndmq_collection,
                            &pattern_list,
                            tag_validated,
                            logging_active,
                            &mut logger,
                            in_records_filename,
                            &mut pattern_hit_counts,
                            &mut nb_hits_tot,
                            &mut nb_records_hit,
                            &mut nb_records_tot,
                            &mut nb_bases,
                            args.filter_matching,
                            args.suppress_output,
                            args.invert_match,
                            &mut writer,
                            &mut json_logger,
                        )?;
                    }
                    Ok(false) => break,
                    Err(e) => anyhow::bail!("Error during BAM record parsing: {}", e),
                }
            }
        }
        "sam" => {
            // Open SAM file for reading
            let mut reader = bam::SamReader::from_path(&args.in_file)
                .with_context(|| format!("Error reading SAM file: {:?}", &args.in_file))?;
            // Get header from SAM file and add program information
            let command_line = env::args().collect::<Vec<String>>().join(" ");
            let mut program_header_line = format!("@PG\tID:{0}\tPN:{0}\tCL:", crate_name!());
            program_header_line.push_str(&command_line);
            program_header_line.push_str(format!("\tVN:{}", crate_version!()).as_str());
            let mut header = reader.header().clone();
            header.push_line(&program_header_line).unwrap();
            // Use empty header if suppress_output is set
            if args.suppress_output {
                header = bam::Header::new();
            }
            // Open file for writing with inferred writer
            let mut writer = match out_file_extension {
                "bam" | "sam" | "STDOUT" => {
                    infer_record_writer(args.threads, &out_file, out_file_extension, header)
                }
                _ => anyhow::bail!("Output file must be a BAM or SAM file."),
            }
            .with_context(|| "Could not create writer.")?;

            // Iterate over SAM records and process them
            let mut record = bam::Record::new();
            loop {
                match reader.read_into(&mut record) {
                    Ok(true) => {
                        process_record(
                            &mut record,
                            ac.as_ref(),
                            &bndmq_collection,
                            &pattern_list,
                            tag_validated,
                            logging_active,
                            &mut logger,
                            in_records_filename,
                            &mut pattern_hit_counts,
                            &mut nb_hits_tot,
                            &mut nb_records_hit,
                            &mut nb_records_tot,
                            &mut nb_bases,
                            args.filter_matching,
                            args.suppress_output,
                            args.invert_match,
                            &mut writer,
                            &mut json_logger,
                        )?;
                    }
                    Ok(false) => break,
                    Err(e) => anyhow::bail!("Error during SAM record parsing: {}", e),
                }
            }
        }
        _ => anyhow::bail!("Input file must be a BAM or SAM file."),
    }

    // Log summary statistics
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
        logger.write_header(&format!("#Total number of hits: {nb_hits_tot}\n"));
        logger.write_header(&format!(
            "#Number of distinct records with a hit: {nb_records_hit}\n"
        ));
        logger.flush();
    }

    // Finalize JSON log if active
    if let Some(jl) = json_logger {
        let input_files_json = serde_json::json!({
            "kmer_file": args.kmer_file.as_ref().map(|p| p.to_string_lossy().to_string()),
            "record_file_1": in_records_filename,
        });
        let pattern_hit_counts_map: HashMap<String, u32> = pattern_list
            .iter()
            .cloned()
            .zip(pattern_hit_counts.iter().copied())
            .collect();
        let meta_information = serde_json::json!({
            "program": crate_name!(),
            "version": crate_version!(),
            "timestamp": Zoned::now().round(Unit::Second)?,
            "subcommand": "tag",
            "command_line": env::args().collect::<Vec<String>>(),
            "search_algorithm": if args.aho_corasick { "Aho-Corasick" } else { "BNDMq" },
            "inverted_matching": args.invert_match,
            "case_insensitive": args.case_insensitive,
            "input_files": input_files_json,
            "tag": from_utf8(&tag_validated).unwrap(),
        });
        let summary_statistics = serde_json::json!({
            "number_of_patterns_searched": pattern_list.len(),
            "number_of_patterns_found": pattern_hit_counts.iter().filter(|&&count| count > 0).count(),
            "number_of_records_searched": nb_records_tot,
            "number_of_characters_searched": nb_bases,
            "number_of_matches": nb_hits_tot,
            "number_of_distinct_records_with_a_hit": nb_records_hit,
        });
        jl.finalize(
            &meta_information,
            &serde_json::json!(pattern_hit_counts_map),
            &summary_statistics,
            None,
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

    // Simple test which does not validate the output!
    #[test]
    fn test_tag_records_no_panics() {
        let in_file = PathBuf::from("tests/data/sample.sam");
        let kmer_seq = vec!["ACGT".to_string()];
        let kmer_file = None;
        let reverse_complement = false;
        let tag = "km".to_string();
        let keep_matching = false;
        let out_log = None;
        let threads = 1;
        let out_file = Some(PathBuf::from("tests/data/sample_tagged.sam"));

        let args = CmdTag {
            in_file,
            kmer_seq: Some(kmer_seq),
            kmer_file,
            reverse_complement,
            canonical: false,
            tag,
            filter_matching: keep_matching,
            out_log,
            json_log: None,
            threads,
            out_file,
            suppress_output: false,
            invert_match: false,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args).unwrap();
        let out_file = PathBuf::from("tests/data/sample_tagged.sam");
        fs::remove_file(out_file).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_tag_records_zero_threads() {
        let in_file = PathBuf::from("tests/data/sample.sam");
        let kmer_seq = vec!["ACGT".to_string()];
        let kmer_file = None;
        let reverse_complement = false;
        let tag = "km".to_string();
        let keep_matching = false;
        let out_log = None;
        let threads = 0;
        let out_file = Some(PathBuf::from("tests/data/sample_tagged.sam"));

        let args = CmdTag {
            in_file,
            kmer_seq: Some(kmer_seq),
            kmer_file,
            reverse_complement,
            canonical: false,
            tag,
            filter_matching: keep_matching,
            out_log,
            json_log: None,
            threads,
            out_file,
            suppress_output: false,
            invert_match: false,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args).unwrap();
        let out_file = PathBuf::from("tests/data/tagged.sam");
        fs::remove_file(out_file).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_tag_records_invalid_tag() {
        let in_file = PathBuf::from("tests/data/sample.sam");
        let kmer_seq = vec!["ACGT".to_string()];
        let kmer_file = None;
        let reverse_complement = false;
        let tag = "kmer".to_string();
        let keep_matching = false;
        let out_log = None;
        let threads = 1;
        let out_file = Some(PathBuf::from("tests/data/tagged.sam"));

        let args = CmdTag {
            in_file,
            kmer_seq: Some(kmer_seq),
            kmer_file,
            reverse_complement,
            canonical: false,
            tag,
            filter_matching: keep_matching,
            out_log,
            json_log: None,
            threads,
            out_file,
            suppress_output: false,
            invert_match: false,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args).unwrap();
        let out_file = PathBuf::from("tests/data/tagged.bam");
        fs::remove_file(out_file).unwrap();
    }

    /// Compare SAM/BAM output with fixture
    fn compare_sam_output(actual_path: &Path, expected_path: &str) -> Result<()> {
        let expected = fs::read_to_string(expected_path)?;
        let actual = fs::read_to_string(actual_path)?;

        // Split into lines
        let expected_lines: Vec<&str> = expected.lines().collect();
        let actual_lines: Vec<&str> = actual.lines().collect();

        // Find the first non-header line in both files
        let expected_record_start = expected_lines
            .iter()
            .position(|line| !line.starts_with('@'))
            .unwrap_or(expected_lines.len());
        let actual_record_start = actual_lines
            .iter()
            .position(|line| !line.starts_with('@'))
            .unwrap_or(actual_lines.len());

        // Compare header lines (excluding @PG)
        let expected_headers: Vec<&str> = expected_lines[..expected_record_start]
            .iter()
            .filter(|line| !line.starts_with("@PG"))
            .copied()
            .collect();
        let actual_headers: Vec<&str> = actual_lines[..actual_record_start]
            .iter()
            .filter(|line| !line.starts_with("@PG"))
            .copied()
            .collect();

        assert_eq!(
            expected_headers, actual_headers,
            "SAM headers (excluding @PG) do not match"
        );

        // Compare record lines
        let expected_records = &expected_lines[expected_record_start..];
        let actual_records = &actual_lines[actual_record_start..];

        assert_eq!(expected_records, actual_records, "SAM records do not match");

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
        assert_eq!(
            expected_json["meta_information"]["tag"], actual_json["meta_information"]["tag"],
            "JSON tag mismatch"
        );

        Ok(())
    }

    // Compare with simple SAM file, including reverse complement and filtering
    // Corresponds to: cargo run -- tag -i tests/fixtures/input/simple.sam -o tests/fixtures/tag/simple.extracted.sam -s CTC -r -l tests/fixtures/tag/simple.log -j tests/fixtures/tag/simple.json -p 2 -m
    #[test]
    fn test_tag_against_sam_fixtures() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_sam = temp_dir.path().join("out.sam");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the tag command
        let args = CmdTag {
            in_file: PathBuf::from("tests/fixtures/input/simple.sam"),
            out_file: Some(out_sam.clone()),
            kmer_seq: Some(vec!["CTC".to_string()]),
            kmer_file: None,
            reverse_complement: true,
            canonical: false,
            tag: "km".to_string(),
            filter_matching: true,
            out_log: Some(out_log.clone()),
            json_log: Some(out_json.clone()),
            threads: 2,
            suppress_output: false,
            invert_match: false,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args)?;

        // Compare outputs with fixtures
        compare_sam_output(&out_sam, "tests/fixtures/tag/simple.extracted.sam")?;
        compare_log_output(&out_log, "tests/fixtures/tag/simple.log")?;
        compare_json_output(&out_json, "tests/fixtures/tag/simple.json")?;

        Ok(())
    }

    // Test inverted matching with simple SAM file
    // Corresponds to: cargo run -- tag -i tests/fixtures/input/simple.sam -o tests/fixtures/tag/simple-inv.extracted.sam -s CTC -r -l tests/fixtures/tag/simple-inv.log -j tests/fixtures/tag/simple-inv.json -p 2 -v
    #[test]
    fn test_tag_against_sam_fixtures_inverted() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_sam = temp_dir.path().join("out.sam");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the tag command with inverted matching
        let args = CmdTag {
            in_file: PathBuf::from("tests/fixtures/input/simple.sam"),
            out_file: Some(out_sam.clone()),
            kmer_seq: Some(vec!["CTC".to_string()]),
            kmer_file: None,
            reverse_complement: true,
            canonical: false,
            tag: "km".to_string(),
            filter_matching: false,
            out_log: Some(out_log.clone()),
            json_log: Some(out_json.clone()),
            threads: 2,
            suppress_output: false,
            invert_match: true,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args)?;

        // Compare outputs with fixtures
        compare_sam_output(&out_sam, "tests/fixtures/tag/simple-inv.extracted.sam")?;
        compare_log_output(&out_log, "tests/fixtures/tag/simple-inv.log")?;
        compare_json_output(&out_json, "tests/fixtures/tag/simple-inv.json")?;

        Ok(())
    }

    // Compare with simple BAM file, including reverse complement but no filtering
    // Corresponds to: cargo run -- tag -i tests/fixtures/input/simple.bam -o tests/fixtures/tag/simple.tagged.extracted.sam -s CTC -r -l tests/fixtures/tag/simple-bam.log -j tests/fixtures/tag/simple-bam.json -p 2
    #[test]
    fn test_tag_against_bam_fixtures() -> Result<()> {
        // Create temporary output files
        let temp_dir = tempfile::tempdir()?;
        let out_sam = temp_dir.path().join("out.sam");
        let out_log = temp_dir.path().join("out.log");
        let out_json = temp_dir.path().join("out.json");

        // Run the tag command
        let args = CmdTag {
            in_file: PathBuf::from("tests/fixtures/input/simple.bam"),
            out_file: Some(out_sam.clone()),
            kmer_seq: Some(vec!["CTC".to_string()]),
            kmer_file: None,
            reverse_complement: true,
            canonical: false,
            tag: "km".to_string(),
            filter_matching: false, // No -m flag
            out_log: Some(out_log.clone()),
            json_log: Some(out_json.clone()),
            threads: 2,
            suppress_output: false,
            invert_match: false,
            q_size: None,
            aho_corasick: false,
            case_insensitive: false,
            lowercase: false,
            uppercase: false,
        };

        tag_records(args)?;

        // Compare outputs with fixtures
        compare_sam_output(&out_sam, "tests/fixtures/tag/simple.tagged.extracted.sam")?;
        compare_log_output(&out_log, "tests/fixtures/tag/simple-bam.log")?;
        compare_json_output(&out_json, "tests/fixtures/tag/simple-bam.json")?;

        Ok(())
    }

    // TODO: Add tests for BAM output - not as easy because of BAM comparison.
}
