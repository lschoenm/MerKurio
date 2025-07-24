//! # Helper functions for working with files and paths.
//!
//! This module contains functions for working with files and paths, such as
//! adding a suffix to a file name, getting the uncompressed file format, and
//! reading k-mers from a file.

use anyhow::{Context, Result};
use needletail::{Sequence, sequence};
use std::{
    fs,
    path::{Path, PathBuf},
    str,
};

/// Checks if the given path is a directory and returns an error if so.
pub fn error_if_directory(path: &Path, description: &str) -> anyhow::Result<()> {
    if path.is_dir() {
        anyhow::bail!(
            "{} '{}' is a directory, not a file.",
            description,
            path.display()
        );
    }
    Ok(())
}

/// Adds provided suffix to file name of a Path, before the first dot.
/// E.g. "sample.fasta.gz" -> "sample_suffix.fasta.gz"
pub fn add_suffix_to_file_prefix(path: &Path, suffix: &str) -> PathBuf {
    let mut pathbuf = path.to_path_buf();
    let file_name = pathbuf
        .file_name()
        .and_then(|name| name.to_str())
        .expect("Invalid file name");

    let mut parts: Vec<String> = file_name.split('.').map(String::from).collect();
    if let Some(first) = parts.first_mut() {
        *first = format!("{first}{suffix}");
    }

    pathbuf.set_file_name(parts.join("."));
    pathbuf
}

/// Get the file format from a file path, ignoring the compression extension if
/// there is one. Returns an error if the file is a directory or has no
/// extension.
pub fn identify_uncompressed_type(path: &Path) -> Result<String> {
    if path.is_dir() {
        anyhow::bail!("The path points to a directory.");
    }

    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .ok_or_else(|| anyhow::anyhow!("Path has no extension"))?;
    let path_decompressed = path.with_extension("");

    if matches!(ext, "gz" | "bz" | "bz2" | "xz") {
        path_decompressed
            .extension()
            .and_then(|e| e.to_str())
            .map(|s| s.to_string())
            .ok_or_else(|| anyhow::anyhow!("Could not determine uncompressed file type"))
    } else {
        Ok(ext.to_string())
    }
}

/// Read k-mers from file to a list (FASTA or one k-mer per line),
/// convert sequence to lowercase or uppercase if flag is set,
/// with or without reverse complements, or compute and use only canonical k-mers.
/// Reverse complement works for IUPAC codes, everything else passes through.
/// TODO: Add more docs about reverse complement and write own implementation with warnings?
/// Returns error if list is empty. Sorts, removes duplicates and empty patterns.
pub fn parse_pattern_list(
    kmer_file: &Option<PathBuf>,
    kmer_seq: Option<Vec<String>>,
    reverse_complement: bool,
    canonical: bool,
    lowercase: bool,
    uppercase: bool,
) -> Result<Vec<String>> {
    // Prioritize reading from path over provided sequence
    let mut pattern_list = match kmer_file {
        Some(path) => read_kmers_from_file(path)
            .with_context(|| format!("Problem reading k-mers from file: {path:?}"))?,
        None => kmer_seq.ok_or_else(|| anyhow::anyhow!("No k-mer sequence provided."))?,
    };

    // Convert sequences to lowercase or uppercase if flag is set
    if lowercase {
        pattern_list = pattern_list.iter().map(|s| s.to_lowercase()).collect();
    } else if uppercase {
        pattern_list = pattern_list.iter().map(|s| s.to_uppercase()).collect();
    }

    // Add reverse complements of k-mers to the list if flag is set
    if reverse_complement {
        let rev_compl_list: Vec<String> = pattern_list
            .iter()
            .map(|pattern| {
                str::from_utf8(&pattern.as_bytes().reverse_complement())
                    .expect("Invalid UTF-8 in reverse complement k-mer.")
                    .to_string()
            })
            .collect();
        pattern_list.extend(rev_compl_list);
    }

    // Use the canonical forms of k-mers
    if canonical {
        pattern_list = pattern_list
            .iter()
            .map(|s| {
                let seq = s.as_bytes();
                let can = sequence::canonical(seq);
                String::from_utf8(can.to_vec()).expect("Invalid UTF-8 in canonical k-mer")
            })
            .collect();
    }

    // Sort pattern list and remove duplicates and empty patterns
    pattern_list.retain(|x| !x.is_empty());
    pattern_list.sort_unstable();
    pattern_list.dedup();

    if pattern_list.is_empty() {
        anyhow::bail!("No k-mers found in file or provided sequence.");
    }

    Ok(pattern_list)
}

/// Read k-mers from file and return them as a vector of strings.
/// Skips empty lines and trims whitespace from the beginning and end of lines.
/// Also skips lines that start with a '#' or '>' character.
/// Returns error if the file contains no k-mers.
pub fn read_kmers_from_file(path: &Path) -> Result<Vec<String>> {
    if path.is_dir() {
        anyhow::bail!(
            "K-mer file path '{}' is a directory, not a file.",
            path.display()
        );
    }

    let content = fs::read_to_string(path).with_context(|| match path.exists() {
        true => format!("Error reading file: {}", path.display()),
        false => "File not found.".to_string(),
    })?;

    let kmer_list: Vec<String> = content
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#') && !line.starts_with('>'))
        .map(|line| line.trim().to_string())
        .collect();

    if kmer_list.is_empty() {
        anyhow::bail!("No k-mers found in the file.");
    }

    Ok(kmer_list)
}

/// Checks for conflicts when writing logs and regular output.
///
/// - Returns an error if both `-l` and `-j` are provided without an argument
///   (i.e. both to stdout).
/// - Returns an error if either log flag is set to stdout while the regular
///   output is also directed to stdout (no `-o` given) and output is not
///   suppressed.
pub fn check_log_flag_conflict(
    out_log: &Option<PathBuf>,
    json_log: &Option<PathBuf>,
    out_file: &Option<PathBuf>,
    suppress_output: bool,
) -> Result<(), String> {
    if let (Some(l), Some(j)) = (out_log, json_log) {
        let l_is_stdout = l.to_str().map(|s| s == "STDOUT").unwrap_or(false);
        let j_is_stdout = j.to_str().map(|s| s == "STDOUT").unwrap_or(false);
        if l_is_stdout && j_is_stdout {
            return Err("Cannot use both -l/--out-log and -j/--json-log with no arguments (both to stdout). Please specify a file for at least one.".to_string());
        }
    }

    let log_to_stdout = out_log
        .as_ref()
        .map(|p| p.to_str().map(|s| s == "STDOUT").unwrap_or(false))
        .unwrap_or(false)
        || json_log
            .as_ref()
            .map(|p| p.to_str().map(|s| s == "STDOUT").unwrap_or(false))
            .unwrap_or(false);

    if log_to_stdout && out_file.is_none() && !suppress_output {
        return Err("Cannot write log to stdout when normal output is also stdout. Specify an output file with -o or suppress output with -S.".to_string());
    }

    Ok(())
}

/// Returns true if Aho-Corasick search algorithm is recommended based on the number of patterns and their length.
pub fn recommend_aho_corasick(pattern_list: &[String]) -> Result<bool> {
    let num_patterns = pattern_list.len();
    let max_len = pattern_list.iter().map(|x| x.len()).max().unwrap();
    if num_patterns >= 14 || max_len > 64 {
        Ok(true)
    } else {
        Ok(false)
    }
}

//
// ---------------------------------- Tests ----------------------------------
//

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_path_is_directory() {
        let path = Path::new("tests/data/");
        let result = error_if_directory(path, "Test path");
        assert!(result.is_err());
    }

    #[test]
    fn test_add_suffix_to_stem_simple() {
        let path = Path::new("tests/data/sample.fasta");
        let new_path = add_suffix_to_file_prefix(path, "_suffix");
        assert_eq!(new_path.to_str().unwrap(), "tests/data/sample_suffix.fasta");
    }

    #[test]
    fn test_add_suffix_to_stem_punctuated() {
        let path = Path::new("tests/data/sample.sorted.fasta");
        let new_path = add_suffix_to_file_prefix(path, "_suffix");
        assert_eq!(
            new_path.to_str().unwrap(),
            "tests/data/sample_suffix.sorted.fasta"
        );
    }

    #[test]
    fn test_get_file_format_simple() {
        let path = Path::new("tests/data/sample.fasta");
        let format = identify_uncompressed_type(path).unwrap();
        assert_eq!(format, "fasta");
    }

    #[test]
    fn test_get_file_format_compressed() {
        let path = Path::new("tests/data/sample.fasta.gz");
        let format = identify_uncompressed_type(path).unwrap();
        assert_eq!(format, "fasta");
    }

    #[test]
    fn test_get_file_format_bz2() {
        let path = Path::new("tests/data/sample.fasta.bz2");
        let format = identify_uncompressed_type(path).unwrap();
        assert_eq!(format, "fasta");
    }

    #[test]
    fn test_get_file_format_xz() {
        let path = Path::new("tests/data/sample.fasta.xz");
        let format = identify_uncompressed_type(path).unwrap();
        assert_eq!(format, "fasta");
    }

    #[test]
    fn test_get_file_format_punctuated() {
        let path = Path::new("tests/data/sample.filtered.fasta.gz");
        let format = identify_uncompressed_type(path).unwrap();
        assert_eq!(format, "fasta");
    }

    #[test]
    fn test_read_kmers_simple() {
        let path = PathBuf::from("tests/data/kmers.txt");
        let kmers = read_kmers_from_file(&path).unwrap();
        assert_eq!(kmers.len(), 3);
        assert!(kmers.contains(&"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string()));
        assert!(kmers.contains(&"AAAATTGCATGAATATTGTAGATCAAAGCACA".to_string()));
        assert!(kmers.contains(&"CTCCGAAGAAGTTGCTGTTCTTGATGGTTATT".to_string()));
    }

    #[test]
    fn test_read_kmers_from_fasta() {
        let path = PathBuf::from("tests/data/kmers.fasta");
        let kmers = read_kmers_from_file(&path).unwrap();
        assert_eq!(kmers.len(), 3);
        assert!(kmers.contains(&"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string()));
        assert!(kmers.contains(&"AAAATTGCATGAATATTGTAGATCAAAGCACA".to_string()));
        assert!(kmers.contains(&"CTCCGAAGAAGTTGCTGTTCTTGATGGTTATT".to_string()));
    }

    #[test]
    fn test_read_kmers_messy() {
        let path = PathBuf::from("tests/data/kmers-messy.txt");
        let kmers = read_kmers_from_file(&path).unwrap();
        assert_eq!(kmers.len(), 3);
        assert!(kmers.contains(&"AAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string()));
        assert!(kmers.contains(&"TTGCATGAATATTGTA".to_string()));
        assert!(kmers.contains(&"CTCCGAAGAAGTTGCTGTTCTTGATGGTTATT".to_string()));
    }

    #[test]
    #[should_panic]
    fn test_read_kmers_empty_file() {
        let path = PathBuf::from("tests/data/kmers-empty.txt");
        read_kmers_from_file(&path).unwrap();
    }

    #[test]
    fn test_parse_pattern_list_from_file_dups() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers-duplicates.txt")),
            None,
            true,
            false,
            false,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list.len(), 4);
    }

    #[test]
    fn test_parse_pattern_list_sorted() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers.txt")),
            None,
            false,
            false,
            false,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list[0], "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(pattern_list[1], "AAAATTGCATGAATATTGTAGATCAAAGCACA");
        assert_eq!(pattern_list[2], "CTCCGAAGAAGTTGCTGTTCTTGATGGTTATT");
    }

    #[test]
    #[should_panic]
    fn test_parse_pattern_list_empty_string() {
        parse_pattern_list(
            &None,
            Some(vec!["".to_string()]),
            false,
            false,
            false,
            false,
        )
        .unwrap();
    }

    #[test]
    fn test_parse_pattern_list_canonical() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers.txt")),
            None,
            false,
            true,
            false,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list.len(), 3);
        assert!(pattern_list.contains(&"AATAACCATCAAGAACAGCAACTTCTTCGGAG".to_string()));
        assert!(pattern_list.contains(&"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string()));
        assert!(pattern_list.contains(&"AAAATTGCATGAATATTGTAGATCAAAGCACA".to_string()));
    }

    #[test]
    fn test_parse_pattern_list_reverse_complement() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers.txt")),
            None,
            true,
            false,
            false,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list.len(), 6);
        assert!(pattern_list.contains(&"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string()));
        assert!(pattern_list.contains(&"AAAATTGCATGAATATTGTAGATCAAAGCACA".to_string()));
        assert!(pattern_list.contains(&"CTCCGAAGAAGTTGCTGTTCTTGATGGTTATT".to_string()));
        assert!(pattern_list.contains(&"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string()));
        assert!(pattern_list.contains(&"TGTGCTTTGATCTACAATATTCATGCAATTTT".to_string()));
        assert!(pattern_list.contains(&"AATAACCATCAAGAACAGCAACTTCTTCGGAG".to_string()));
    }

    #[test]
    fn test_parse_pattern_list_aa() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers-aa.txt")),
            None,
            false,
            false,
            false,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list.len(), 3);
        assert!(pattern_list.contains(&"MDLQENLVSDAGDDHMV".to_string()));
        assert!(pattern_list.contains(&"DIVVEPHSNRDIGIVDE".to_string()));
        assert!(pattern_list.contains(&"FNIGGDVGFSGDLDLEP".to_string()));
    }

    #[test]
    fn test_parse_pattern_list_lowercase() {
        let pattern_list = parse_pattern_list(
            &Some(PathBuf::from("tests/data/kmers.txt")),
            None,
            false,
            false,
            true,
            false,
        )
        .unwrap();
        assert_eq!(pattern_list.len(), 3);
        assert!(pattern_list.contains(&"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa".to_string()));
        assert!(pattern_list.contains(&"aaaattgcatgaatattgtagatcaaagcaca".to_string()));
        assert!(pattern_list.contains(&"ctccgaagaagttgctgttcttgatggttatt".to_string()));
    }

    #[test]
    fn test_logs_both_none() {
        assert!(check_log_flag_conflict(&None, &None, &None, false).is_ok());
    }

    #[test]
    fn test_logs_one_stdout_one_none() {
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("STDOUT")),
                &None,
                &Some(PathBuf::from("file.fa")),
                false
            )
            .is_ok()
        );
        assert!(
            check_log_flag_conflict(
                &None,
                &Some(PathBuf::from("STDOUT")),
                &Some(PathBuf::from("file.fa")),
                false
            )
            .is_ok()
        );
    }

    #[test]
    fn test_logs_one_file_one_none() {
        assert!(
            check_log_flag_conflict(&Some(PathBuf::from("file.log")), &None, &None, false).is_ok()
        );
        assert!(
            check_log_flag_conflict(&None, &Some(PathBuf::from("file.json")), &None, false).is_ok()
        );
    }

    #[test]
    fn test_logs_both_stdout_conflict() {
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("STDOUT")),
                &Some(PathBuf::from("STDOUT")),
                &None,
                false
            )
            .is_err()
        );
    }

    #[test]
    fn test_logs_both_file() {
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("file1.log")),
                &Some(PathBuf::from("file2.json")),
                &None,
                false
            )
            .is_ok()
        );
    }

    #[test]
    fn test_logs_one_stdout_two_files() {
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("STDOUT")),
                &Some(PathBuf::from("file.json")),
                &Some(PathBuf::from("file.fa")),
                false
            )
            .is_ok()
        );
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("file.log")),
                &Some(PathBuf::from("STDOUT")),
                &Some(PathBuf::from("file.fa")),
                false
            )
            .is_ok()
        );
    }

    #[test]
    fn test_logs_one_stdout_one_file() {
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("STDOUT")),
                &Some(PathBuf::from("file.json")),
                &None,
                false
            )
            .is_err()
        );
        assert!(
            check_log_flag_conflict(
                &Some(PathBuf::from("file.log")),
                &Some(PathBuf::from("STDOUT")),
                &None,
                false
            )
            .is_err()
        );
    }

    #[test]
    fn test_logs_conflict_with_stdout_output() {
        assert!(
            check_log_flag_conflict(&Some(PathBuf::from("STDOUT")), &None, &None, false).is_err()
        );
        assert!(
            check_log_flag_conflict(&None, &Some(PathBuf::from("STDOUT")), &None, false).is_err()
        );
        // No conflict when output is suppressed
        assert!(
            check_log_flag_conflict(&Some(PathBuf::from("STDOUT")), &None, &None, true).is_ok()
        );
    }

    #[test]
    fn test_tune_search_algorithm_patterns_few() {
        let ac = recommend_aho_corasick(&["AAA".to_string(), "CCC".to_string()]).unwrap();
        assert_eq!(ac, false);
    }

    #[test]
    fn test_tune_search_algorithm_patterns_long() {
        let ac = recommend_aho_corasick(&[
            "AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTAAAAAAAACCCCCCCCGGGGGGGGTTTTTTTTA".to_string(),
        ])
        .unwrap();
        assert_eq!(ac, true);
    }
}
