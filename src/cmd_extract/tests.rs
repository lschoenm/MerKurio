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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
    let actual_json: serde_json::Value = serde_json::from_str(&fs::read_to_string(actual_path)?)?;

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

fn assert_single_runs_match(actual: &SingleExtractRun, expected: &SingleExtractRun) -> Result<()> {
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
    let actual_json: serde_json::Value = serde_json::from_str(&fs::read_to_string(actual_path)?)?;

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
        expected_json["paired_end_reads_statistics"], actual_json["paired_end_reads_statistics"],
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
    let parallel_2 = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 2))?;
    let parallel_4 = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;

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
    let parallel = run_single_extract(SingleExtractOptions::simple(vec!["ACG".to_string()], 4))?;

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
            effective_total_threads: 8,
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
            effective_total_threads: 8,
            auto_detected: false,
            clamped: true,
        }
    );
}

#[test]
fn test_extract_total_thread_budget_reserves_one_reader_thread() {
    assert_eq!(matching_thread_count(1), 0);
    assert_eq!(matching_thread_count(2), 1);
    assert_eq!(matching_thread_count(4), 3);
}

#[test]
fn test_extract_rejects_zero_chunk_size() {
    let args = CmdExtract {
        in_fastx: PathBuf::from("tests/fixtures/input/simple.fasta"),
        in_fastq_2: None,
        kmer_seq: Some(vec!["ACG".to_string()]),
        kmer_file: None,
        out_fastx: None,
        q_size: None,
        aho_corasick: false,
        reverse_complement: true,
        canonical: false,
        out_log: None,
        suppress_output: false,
        json_log: None,
        invert_match: false,
        case_insensitive: false,
        lowercase: false,
        uppercase: false,
        threads: 4,
        chunk_size: 0,
    };

    let error = extract_records(args).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Extract chunk size must be greater than zero")
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
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
        chunk_size: DEFAULT_EXTRACT_PARALLEL_CHUNK_SIZE,
    };

    let error = extract_records(args).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("The two input files have a different number of records")
    );

    Ok(())
}
