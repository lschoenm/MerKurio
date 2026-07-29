#![allow(dead_code)]
#![allow(clippy::upper_case_acronyms)]

use std::collections::HashSet;
use std::env;
use std::fs::{self, File, OpenOptions};
use std::hint::black_box;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

#[path = "../../../src/pattern_matching.rs"]
mod pattern_matching;
#[path = "../../../src/pattern_preprocessing.rs"]
mod pattern_preprocessing;

use pattern_matching::PatternMatcher;

const DEFAULT_SEED: u64 = 0x4d65_724b_7572_696f;
const DEFAULT_K_VALUES: &[usize] = &[4, 6, 8, 10, 12, 16, 20, 24, 31, 40, 48, 64, 80, 96, 128];
const DEFAULT_PATTERN_COUNTS: &[usize] = &[
    1, 2, 4, 8, 16, 32, 64, 128, 256, 1_000, 10_000, 100_000, 1_000_000,
];

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum Algorithm {
    Bndmq,
    Hash,
    AhoCorasick,
}

impl Algorithm {
    const ALL: [Self; 3] = [Self::Bndmq, Self::Hash, Self::AhoCorasick];

    fn name(self) -> &'static str {
        match self {
            Self::Bndmq => "bndmq",
            Self::Hash => "hash",
            Self::AhoCorasick => "aho_corasick",
        }
    }

    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "bndmq" => Ok(Self::Bndmq),
            "hash" => Ok(Self::Hash),
            "aho_corasick" => Ok(Self::AhoCorasick),
            _ => Err(format!("Unknown algorithm: {value}")),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum SearchMode {
    First,
    All,
}

impl SearchMode {
    const ALL: [Self; 2] = [Self::First, Self::All];

    fn name(self) -> &'static str {
        match self {
            Self::First => "first",
            Self::All => "all",
        }
    }

    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "first" => Ok(Self::First),
            "all" => Ok(Self::All),
            _ => Err(format!("Unknown search mode: {value}")),
        }
    }
}

#[derive(Clone)]
struct Config {
    k_values: Vec<usize>,
    pattern_counts: Vec<usize>,
    algorithms: Vec<Algorithm>,
    modes: Vec<SearchMode>,
    sequence_length: usize,
    sequences_per_bank: usize,
    target_patterns_per_cell: usize,
    max_banks: usize,
    runs: usize,
    target_run_ms: u64,
    max_sample_ms: u64,
    cell_timeout_seconds: u64,
    seed: u64,
    output: PathBuf,
    status_output: PathBuf,
    metadata: PathBuf,
    worker: Option<WorkerConfig>,
}

#[derive(Clone)]
struct WorkerConfig {
    k: usize,
    pattern_count: usize,
    algorithm: Algorithm,
    mode: SearchMode,
    output: PathBuf,
}

impl Default for Config {
    fn default() -> Self {
        let results = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
        Self {
            k_values: DEFAULT_K_VALUES.to_vec(),
            pattern_counts: DEFAULT_PATTERN_COUNTS.to_vec(),
            algorithms: Algorithm::ALL.to_vec(),
            modes: SearchMode::ALL.to_vec(),
            sequence_length: 150,
            sequences_per_bank: 256,
            target_patterns_per_cell: 256,
            max_banks: 64,
            runs: 12,
            target_run_ms: 20,
            max_sample_ms: 1_000,
            cell_timeout_seconds: 30,
            seed: DEFAULT_SEED,
            output: results.join("algorithm_sweep.csv"),
            status_output: results.join("cell_status.csv"),
            metadata: results.join("metadata.txt"),
            worker: None,
        }
    }
}

impl Config {
    fn parse() -> Result<Self, String> {
        let mut config = Self::default();
        let mut worker = false;
        let mut worker_k = None;
        let mut worker_patterns = None;
        let mut worker_algorithm = None;
        let mut worker_mode = None;
        let mut worker_output = None;
        let mut args = env::args().skip(1);

        while let Some(flag) = args.next() {
            let value = |args: &mut std::iter::Skip<std::env::Args>,
                         flag: &str|
             -> Result<String, String> {
                args.next()
                    .ok_or_else(|| format!("Missing value for {flag}."))
            };
            match flag.as_str() {
                "--k-values" => {
                    config.k_values = parse_usize_list(&value(&mut args, &flag)?, &flag)?
                }
                "--pattern-counts" => {
                    config.pattern_counts = parse_usize_list(&value(&mut args, &flag)?, &flag)?
                }
                "--algorithms" => {
                    config.algorithms = value(&mut args, &flag)?
                        .split(',')
                        .map(|item| Algorithm::parse(item.trim()))
                        .collect::<Result<_, _>>()?
                }
                "--modes" => {
                    config.modes = value(&mut args, &flag)?
                        .split(',')
                        .map(|item| SearchMode::parse(item.trim()))
                        .collect::<Result<_, _>>()?
                }
                "--sequence-length" => {
                    config.sequence_length = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--sequences" => {
                    config.sequences_per_bank = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--target-patterns-per-cell" => {
                    config.target_patterns_per_cell = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--max-banks" => config.max_banks = parse_usize(&value(&mut args, &flag)?, &flag)?,
                "--runs" | "--samples" => {
                    config.runs = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--target-ms" => {
                    config.target_run_ms = parse_u64(&value(&mut args, &flag)?, &flag)?
                }
                "--max-sample-ms" => {
                    config.max_sample_ms = parse_u64(&value(&mut args, &flag)?, &flag)?
                }
                "--cell-timeout-seconds" => {
                    config.cell_timeout_seconds = parse_u64(&value(&mut args, &flag)?, &flag)?
                }
                "--seed" => config.seed = parse_u64(&value(&mut args, &flag)?, &flag)?,
                "--output" => config.output = PathBuf::from(value(&mut args, &flag)?),
                "--status-output" => config.status_output = PathBuf::from(value(&mut args, &flag)?),
                "--metadata" => config.metadata = PathBuf::from(value(&mut args, &flag)?),
                "--worker" => worker = true,
                "--worker-k" => worker_k = Some(parse_usize(&value(&mut args, &flag)?, &flag)?),
                "--worker-patterns" => {
                    worker_patterns = Some(parse_usize(&value(&mut args, &flag)?, &flag)?)
                }
                "--worker-algorithm" => {
                    worker_algorithm = Some(Algorithm::parse(&value(&mut args, &flag)?)?)
                }
                "--worker-mode" => {
                    worker_mode = Some(SearchMode::parse(&value(&mut args, &flag)?)?)
                }
                "--worker-output" => worker_output = Some(PathBuf::from(value(&mut args, &flag)?)),
                "-h" | "--help" => {
                    print_help();
                    std::process::exit(0);
                }
                _ => return Err(format!("Unknown option: {flag}. Use --help for usage.")),
            }
        }

        if worker {
            config.worker = Some(WorkerConfig {
                k: worker_k.ok_or("Worker is missing --worker-k.")?,
                pattern_count: worker_patterns.ok_or("Worker is missing --worker-patterns.")?,
                algorithm: worker_algorithm.ok_or("Worker is missing --worker-algorithm.")?,
                mode: worker_mode.ok_or("Worker is missing --worker-mode.")?,
                output: worker_output.ok_or("Worker is missing --worker-output.")?,
            });
        }
        config.validate()?;
        Ok(config)
    }

    fn validate(&mut self) -> Result<(), String> {
        self.k_values.sort_unstable();
        self.k_values.dedup();
        self.pattern_counts.sort_unstable();
        self.pattern_counts.dedup();

        deduplicate(&mut self.algorithms);
        deduplicate(&mut self.modes);
        if self.k_values.is_empty()
            || self.pattern_counts.is_empty()
            || self.algorithms.is_empty()
            || self.modes.is_empty()
        {
            return Err(
                "k-values, pattern-counts, algorithms, and modes cannot be empty.".to_string(),
            );
        }
        if *self.pattern_counts.last().unwrap() > 1_000_000 {
            return Err("pattern-counts currently supports at most 1000000 patterns.".to_string());
        }
        if self.k_values[0] == 0 {
            return Err("Pattern lengths must be positive.".to_string());
        }
        if *self.k_values.last().unwrap() > self.sequence_length {
            return Err("Pattern lengths cannot exceed sequence-length.".to_string());
        }
        if self.pattern_counts[0] == 0
            || self.sequence_length == 0
            || self.sequences_per_bank == 0
            || self.target_patterns_per_cell == 0
            || self.max_banks == 0
            || self.runs == 0
            || self.target_run_ms == 0
            || self.max_sample_ms == 0
            || self.cell_timeout_seconds == 0
        {
            return Err("All counts and timing limits must be positive.".to_string());
        }
        if !self.sequences_per_bank.is_multiple_of(2) {
            return Err("sequences must be even for balanced insertion.".to_string());
        }
        Ok(())
    }

    fn bank_count(&self, pattern_count: usize) -> usize {
        self.target_patterns_per_cell
            .div_ceil(pattern_count)
            .clamp(1, self.max_banks)
    }
}

fn parse_usize(value: &str, flag: &str) -> Result<usize, String> {
    value
        .parse()
        .map_err(|_| format!("Invalid integer for {flag}: {value}"))
}

fn parse_u64(value: &str, flag: &str) -> Result<u64, String> {
    value
        .parse()
        .map_err(|_| format!("Invalid integer for {flag}: {value}"))
}

fn parse_usize_list(value: &str, flag: &str) -> Result<Vec<usize>, String> {
    value
        .split(',')
        .map(|item| parse_usize(item.trim(), flag))
        .collect()
}

fn deduplicate<T: Copy + Eq + std::hash::Hash>(values: &mut Vec<T>) {
    let mut seen = HashSet::new();
    values.retain(|value| seen.insert(*value));
}

fn print_help() {
    println!(
        "\
MerKurio pattern-matching algorithm selection tuning

Options:
  --k-values LIST              Comma-separated pattern lengths
  --pattern-counts LIST        Comma-separated pattern counts, up to 1000000
  --algorithms LIST            bndmq,hash,aho_corasick [default: all]
  --modes LIST                 first,all [default: both]
  --sequence-length N          Bases per sequence [default: 150]
  --sequences N                Sequences per independent bank [default: 256]
  --target-patterns-per-cell N Target aggregate pattern count at small p [default: 256]
  --max-banks N                Maximum independent pattern banks [default: 64]
  --runs N                     Measurement rounds [default: 12]
  --target-ms N                Approximate minimum duration per timed run [default: 20]
  --max-sample-ms N            Reject a cell if one search iteration exceeds this [default: 1000]
  --cell-timeout-seconds N     Hard timeout for a worker cell [default: 30]
  --seed N                     Deterministic unsigned integer seed
  --output PATH                Raw successful timing CSV
  --status-output PATH         Cell status CSV
  --metadata PATH              Run metadata path
  -h, --help                   Show this help"
    );
}

#[derive(Clone)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    fn range(&mut self, upper: usize) -> usize {
        (self.next_u64() % upper as u64) as usize
    }

    fn dna_base(&mut self) -> u8 {
        b"ACGT"[self.range(4)]
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            values.swap(index, self.range(index + 1));
        }
    }
}

fn random_dna(rng: &mut SplitMix64, length: usize) -> Vec<u8> {
    (0..length).map(|_| rng.dna_base()).collect()
}

fn maximum_distinct_patterns(k: usize) -> usize {
    4usize.checked_pow(k as u32).unwrap_or(usize::MAX)
}

fn uniqueness_digits(pattern_count: usize) -> usize {
    let mut digits = 0;
    let mut capacity = 1usize;
    while capacity < pattern_count {
        capacity = capacity.saturating_mul(4);
        digits += 1;
    }
    digits.max(1)
}

fn build_patterns(seed: u64, k: usize, pattern_count: usize) -> Vec<String> {
    let digits = uniqueness_digits(pattern_count);
    let modulus = 1usize << (digits * 2);
    let multiplier = 0x9e37_79b1usize | 1;
    let offset = (seed as usize) & (modulus - 1);
    let mut patterns = Vec::with_capacity(pattern_count);

    for pattern_index in 0..pattern_count {
        let mut rng =
            SplitMix64::new(seed ^ (pattern_index as u64).wrapping_mul(0xd6e8_feb8_6659_fd93));
        let mut pattern = random_dna(&mut rng, k);
        let mut code = pattern_index.wrapping_mul(multiplier).wrapping_add(offset) & (modulus - 1);
        for digit in 0..digits {
            let position = digit * k / digits;
            pattern[position] = b"ACGT"[code & 3];
            code >>= 2;
        }
        patterns.push(String::from_utf8(pattern).unwrap());
    }
    patterns
}

fn build_matcher(patterns: &[String], algorithm: Algorithm) -> Result<PatternMatcher, String> {
    let result = match algorithm {
        Algorithm::Bndmq => PatternMatcher::new(patterns, false, false, false, None),
        Algorithm::Hash => PatternMatcher::new(patterns, false, true, false, None),
        Algorithm::AhoCorasick => PatternMatcher::new(patterns, true, false, false, None),
    };
    result.map_err(|error| error.to_string())
}

struct MatcherBank {
    matcher: PatternMatcher,
    sequences: Vec<Vec<u8>>,
}

fn build_banks(config: &Config, worker: &WorkerConfig) -> Result<Vec<MatcherBank>, String> {
    let bank_count = config.bank_count(worker.pattern_count);
    let mut banks = Vec::with_capacity(bank_count);

    for bank_index in 0..bank_count {
        let bank_seed = config.seed
            ^ (worker.k as u64).wrapping_mul(0xa076_1d64_78bd_642f)
            ^ (worker.pattern_count as u64).wrapping_mul(0xe703_7ed1_a0b4_28db)
            ^ (bank_index as u64).wrapping_mul(0x8ebc_6af0_9c88_c6e3);
        let patterns = build_patterns(bank_seed, worker.k, worker.pattern_count);
        let matcher = build_matcher(&patterns, worker.algorithm)?;
        let mut rng = SplitMix64::new(bank_seed ^ 0x5899_65cc_7537_4cc3);
        let mut sequences = Vec::with_capacity(config.sequences_per_bank);
        for sequence_index in 0..config.sequences_per_bank {
            let mut sequence = random_dna(&mut rng, config.sequence_length);
            if sequence_index.is_multiple_of(2) {
                let pattern = patterns[rng.range(patterns.len())].as_bytes();
                let position = rng.range(config.sequence_length - worker.k + 1);
                sequence[position..position + worker.k].copy_from_slice(pattern);
            }
            sequences.push(sequence);
        }
        banks.push(MatcherBank { matcher, sequences });
    }
    Ok(banks)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SearchResult {
    checksum: u64,
    matching_records: u64,
    matches: u64,
}

fn mix_event(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn search_banks(
    banks: &[MatcherBank],
    bank_order: &[usize],
    mode: SearchMode,
    iterations: usize,
) -> SearchResult {
    let mut result = SearchResult {
        checksum: 0,
        matching_records: 0,
        matches: 0,
    };
    for _ in 0..iterations {
        for &bank_index in bank_order {
            let bank = &banks[bank_index];
            for (sequence_index, sequence) in bank.sequences.iter().enumerate() {
                let sequence = black_box(sequence.as_slice());
                match mode {
                    SearchMode::First => {
                        if bank.matcher.find_any(sequence) {
                            result.matching_records += 1;
                            result.matches += 1;
                            result.checksum = result.checksum.wrapping_add(mix_event(
                                ((bank_index as u64) << 32) ^ sequence_index as u64,
                            ));
                        }
                    }
                    SearchMode::All => {
                        let mut matched = false;
                        bank.matcher.for_each_match(sequence, |hit| {
                            matched = true;
                            result.matches += 1;
                            result.checksum = result.checksum.wrapping_add(mix_event(
                                ((bank_index as u64) << 48)
                                    ^ ((sequence_index as u64) << 32)
                                    ^ ((hit.pattern_index as u64) << 16)
                                    ^ hit.position as u64,
                            ));
                        });
                        if matched {
                            result.matching_records += 1;
                        }
                    }
                }
            }
        }
    }
    black_box(result)
}

fn timed_search(
    banks: &[MatcherBank],
    bank_order: &[usize],
    mode: SearchMode,
    iterations: usize,
) -> (u128, SearchResult) {
    let start = Instant::now();
    let result = search_banks(banks, bank_order, mode, iterations);
    (start.elapsed().as_nanos(), result)
}

enum WorkerOutcome {
    Complete,
    SoftTimeout(String),
}

fn run_worker(config: &Config, worker: &WorkerConfig) -> Result<WorkerOutcome, String> {
    let build_start = Instant::now();
    let banks = build_banks(config, worker)?;
    let build_ns = build_start.elapsed().as_nanos();
    let natural_order: Vec<usize> = (0..banks.len()).collect();
    let validation = search_banks(&banks, &natural_order, worker.mode, 1);
    let (calibration_ns, _) = timed_search(&banks, &natural_order, worker.mode, 1);
    let calibration_ms = calibration_ns as f64 / 1_000_000.0;
    if calibration_ns > u128::from(config.max_sample_ms) * 1_000_000 {
        return Ok(WorkerOutcome::SoftTimeout(format!(
            "single iteration took {calibration_ms:.3} ms"
        )));
    }

    let target_ns = u128::from(config.target_run_ms) * 1_000_000;
    let iterations = target_ns.div_ceil(calibration_ns.max(1)).clamp(1, 100_000) as usize;
    black_box(search_banks(
        &banks,
        &natural_order,
        worker.mode,
        iterations,
    ));

    if let Some(parent) = worker.output.parent() {
        fs::create_dir_all(parent).map_err(|error| error.to_string())?;
    }
    let mut output =
        BufWriter::new(File::create(&worker.output).map_err(|error| error.to_string())?);
    let bases_per_iteration = banks
        .iter()
        .map(|bank| bank.sequences.iter().map(Vec::len).sum::<usize>())
        .sum::<usize>();
    let mut round_rng = SplitMix64::new(
        config.seed
            ^ (worker.k as u64).wrapping_mul(0x2d35_8dcc_aa6c_78a5)
            ^ (worker.pattern_count as u64).wrapping_mul(0x8bb8_4b93_962e_acc9)
            ^ match worker.mode {
                SearchMode::First => 0x4f74_4522_5e82_c488,
                SearchMode::All => 0x9e6c_63d0_676a_9a99,
            },
    );

    for run_index in 0..config.runs {
        let mut bank_order = natural_order.clone();
        round_rng.shuffle(&mut bank_order);
        let (elapsed_ns, timed_result) = timed_search(&banks, &bank_order, worker.mode, iterations);
        let bases = bases_per_iteration.saturating_mul(iterations);
        let ns_per_base = elapsed_ns as f64 / bases as f64;
        writeln!(
            output,
            "{},{},{},{},{run_index},{iterations},{},{},{},{build_ns},{bases},{elapsed_ns},{ns_per_base:.9},{},{},{},{}",
            worker.k,
            worker.pattern_count,
            worker.algorithm.name(),
            worker.mode.name(),
            banks.len(),
            config.sequences_per_bank,
            config.sequence_length,
            validation.checksum,
            validation.matching_records,
            validation.matches,
            timed_result.checksum,
        )
        .map_err(|error| error.to_string())?;
    }
    output.flush().map_err(|error| error.to_string())?;
    Ok(WorkerOutcome::Complete)
}

fn command_output(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|output| output.status.success())
        .map(|output| String::from_utf8_lossy(&output.stdout).trim().to_string())
        .unwrap_or_else(|| "unavailable".to_string())
}

fn join_usize(values: &[usize]) -> String {
    values
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn write_metadata(config: &Config) -> std::io::Result<()> {
    if let Some(parent) = config.metadata.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut output = BufWriter::new(File::create(&config.metadata)?);
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    writeln!(output, "unix_timestamp={timestamp}")?;
    writeln!(
        output,
        "git_commit={}",
        command_output("git", &["rev-parse", "HEAD"])
    )?;
    writeln!(
        output,
        "git_status={}",
        command_output("git", &["status", "--short", "--untracked-files=all"]).replace('\n', "; ")
    )?;
    writeln!(
        output,
        "rustc={}",
        command_output("rustc", &["--version", "--verbose"]).replace('\n', "; ")
    )?;
    writeln!(output, "system={}", command_output("uname", &["-a"]))?;
    writeln!(output, "corpus=mixed_inserted_50_percent")?;
    writeln!(output, "k_values={}", join_usize(&config.k_values))?;
    writeln!(
        output,
        "pattern_counts={}",
        join_usize(&config.pattern_counts)
    )?;
    writeln!(
        output,
        "algorithms={}",
        config
            .algorithms
            .iter()
            .map(|algorithm| algorithm.name())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    writeln!(
        output,
        "modes={}",
        config
            .modes
            .iter()
            .map(|mode| mode.name())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    writeln!(output, "sequence_length={}", config.sequence_length)?;
    writeln!(output, "sequences_per_bank={}", config.sequences_per_bank)?;
    writeln!(
        output,
        "target_patterns_per_cell={}",
        config.target_patterns_per_cell
    )?;
    writeln!(output, "max_banks={}", config.max_banks)?;
    writeln!(output, "runs={}", config.runs)?;
    writeln!(output, "target_run_ms={}", config.target_run_ms)?;
    writeln!(output, "max_sample_ms={}", config.max_sample_ms)?;
    writeln!(
        output,
        "cell_timeout_seconds={}",
        config.cell_timeout_seconds
    )?;
    writeln!(output, "seed={}", config.seed)?;
    Ok(())
}

fn csv_text(value: &str) -> String {
    format!("\"{}\"", value.replace('"', "\"\"").replace('\n', " "))
}

fn append_file(source: &Path, destination: &Path) -> std::io::Result<()> {
    let bytes = fs::read(source)?;
    let mut output = OpenOptions::new().append(true).open(destination)?;
    output.write_all(&bytes)
}

struct CellResult {
    status: &'static str,
    wall_ms: u128,
    detail: String,
}

fn run_worker_process(
    config: &Config,
    worker: &WorkerConfig,
    temporary_output: &Path,
    temporary_stderr: &Path,
) -> Result<CellResult, String> {
    let executable = env::current_exe().map_err(|error| error.to_string())?;
    let stderr = File::create(temporary_stderr).map_err(|error| error.to_string())?;
    let mut child = Command::new(executable)
        .arg("--worker")
        .args(["--worker-k", &worker.k.to_string()])
        .args(["--worker-patterns", &worker.pattern_count.to_string()])
        .args(["--worker-algorithm", worker.algorithm.name()])
        .args(["--worker-mode", worker.mode.name()])
        .args(["--worker-output", &temporary_output.display().to_string()])
        .args(["--sequence-length", &config.sequence_length.to_string()])
        .args(["--sequences", &config.sequences_per_bank.to_string()])
        .args([
            "--target-patterns-per-cell",
            &config.target_patterns_per_cell.to_string(),
        ])
        .args(["--max-banks", &config.max_banks.to_string()])
        .args(["--runs", &config.runs.to_string()])
        .args(["--target-ms", &config.target_run_ms.to_string()])
        .args(["--max-sample-ms", &config.max_sample_ms.to_string()])
        .args([
            "--cell-timeout-seconds",
            &config.cell_timeout_seconds.to_string(),
        ])
        .args(["--seed", &config.seed.to_string()])
        .stdout(Stdio::null())
        .stderr(Stdio::from(stderr))
        .spawn()
        .map_err(|error| error.to_string())?;

    let start = Instant::now();
    let timeout = Duration::from_secs(config.cell_timeout_seconds);
    let exit_status = loop {
        if let Some(status) = child.try_wait().map_err(|error| error.to_string())? {
            break Some(status);
        }
        if start.elapsed() >= timeout {
            child.kill().map_err(|error| error.to_string())?;
            child.wait().map_err(|error| error.to_string())?;
            break None;
        }
        thread::sleep(Duration::from_millis(25));
    };
    let wall_ms = start.elapsed().as_millis();
    let stderr = fs::read_to_string(temporary_stderr).unwrap_or_default();

    match exit_status {
        None => Ok(CellResult {
            status: "timed_out",
            wall_ms,
            detail: format!("hard timeout after {} seconds", config.cell_timeout_seconds),
        }),
        Some(status) if status.success() => Ok(CellResult {
            status: "ok",
            wall_ms,
            detail: String::new(),
        }),
        Some(status) if status.code() == Some(3) => Ok(CellResult {
            status: "timed_out",
            wall_ms,
            detail: stderr.trim().to_string(),
        }),
        Some(status) => Ok(CellResult {
            status: "error",
            wall_ms,
            detail: format!("worker exit {status}: {}", stderr.trim()),
        }),
    }
}

fn write_status(
    output: &mut impl Write,
    k: usize,
    patterns: usize,
    algorithm: Algorithm,
    mode: SearchMode,
    result: &CellResult,
) -> std::io::Result<()> {
    writeln!(
        output,
        "{k},{patterns},{},{},{},{},{}",
        algorithm.name(),
        mode.name(),
        result.status,
        result.wall_ms,
        csv_text(&result.detail),
    )
}

fn run_parent(config: &Config) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(parent) = config.output.parent() {
        fs::create_dir_all(parent)?;
    }
    write_metadata(config)?;
    let mut raw = BufWriter::new(File::create(&config.output)?);
    writeln!(
        raw,
        "k,patterns,algorithm,mode,run,iterations,banks,sequences_per_bank,sequence_length,build_ns,bases,elapsed_ns,ns_per_base,validation_checksum,matching_records,matches,timed_checksum"
    )?;
    raw.flush()?;
    let mut statuses = BufWriter::new(File::create(&config.status_output)?);
    writeln!(statuses, "k,patterns,algorithm,mode,status,wall_ms,detail")?;

    let temporary_directory = config
        .output
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(".workers");
    fs::create_dir_all(&temporary_directory)?;
    let mut pruned = HashSet::new();

    for &k in &config.k_values {
        for &mode in &config.modes {
            for &pattern_count in &config.pattern_counts {
                let mut algorithms = config.algorithms.clone();
                let mut order_rng = SplitMix64::new(
                    config.seed
                        ^ (k as u64).wrapping_mul(0x243f_6a88_85a3_08d3)
                        ^ (pattern_count as u64).wrapping_mul(0x1319_8a2e_0370_7344)
                        ^ match mode {
                            SearchMode::First => 0xa409_3822_299f_31d0,
                            SearchMode::All => 0x082e_fa98_ec4e_6c89,
                        },
                );
                order_rng.shuffle(&mut algorithms);

                for algorithm in algorithms {
                    let worker = WorkerConfig {
                        k,
                        pattern_count,
                        algorithm,
                        mode,
                        output: PathBuf::new(),
                    };
                    let key = (k, mode, algorithm);
                    let invalid_detail = if pattern_count > maximum_distinct_patterns(k) {
                        Some(format!(
                            "only {} distinct DNA patterns exist",
                            maximum_distinct_patterns(k)
                        ))
                    } else if algorithm == Algorithm::Bndmq && k > usize::BITS as usize {
                        Some(format!("BNDMq limit is {} bytes", usize::BITS))
                    } else {
                        None
                    };

                    if let Some(detail) = invalid_detail {
                        let result = CellResult {
                            status: "invalid",
                            wall_ms: 0,
                            detail,
                        };
                        write_status(&mut statuses, k, pattern_count, algorithm, mode, &result)?;
                        continue;
                    }
                    if pruned.contains(&key) {
                        let result = CellResult {
                            status: "pruned",
                            wall_ms: 0,
                            detail: "a smaller pattern count timed out".to_string(),
                        };
                        write_status(&mut statuses, k, pattern_count, algorithm, mode, &result)?;
                        continue;
                    }

                    eprintln!(
                        "k={k} patterns={pattern_count} mode={} algorithm={}",
                        mode.name(),
                        algorithm.name()
                    );
                    let stem = format!("{k}-{pattern_count}-{}-{}", algorithm.name(), mode.name());
                    let worker_output = temporary_directory.join(format!("{stem}.csv"));
                    let worker_stderr = temporary_directory.join(format!("{stem}.stderr"));
                    let result =
                        run_worker_process(config, &worker, &worker_output, &worker_stderr)?;
                    if result.status == "ok" {
                        append_file(&worker_output, &config.output)?;
                    } else if matches!(result.status, "timed_out" | "error") {
                        pruned.insert(key);
                    }
                    write_status(&mut statuses, k, pattern_count, algorithm, mode, &result)?;
                    statuses.flush()?;
                    let _ = fs::remove_file(worker_output);
                    let _ = fs::remove_file(worker_stderr);
                }
            }
        }
    }
    Ok(())
}

fn main() {
    let config = Config::parse().unwrap_or_else(|error| {
        eprintln!("Error: {error}");
        std::process::exit(2);
    });

    if let Some(worker) = &config.worker {
        match run_worker(&config, worker) {
            Ok(WorkerOutcome::Complete) => {}
            Ok(WorkerOutcome::SoftTimeout(detail)) => {
                eprintln!("{detail}");
                std::process::exit(3);
            }
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        }
    } else if let Err(error) = run_parent(&config) {
        eprintln!("Tuning failed: {error:#}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tuning_tests {
    use super::*;

    fn small_config() -> Config {
        Config {
            sequence_length: 40,
            sequences_per_bank: 8,
            target_patterns_per_cell: 8,
            max_banks: 2,
            runs: 1,
            target_run_ms: 1,
            max_sample_ms: 100,
            ..Config::default()
        }
    }

    #[test]
    fn generated_patterns_are_deterministic_and_distinct() {
        let first = build_patterns(DEFAULT_SEED, 10, 1_000);
        let second = build_patterns(DEFAULT_SEED, 10, 1_000);
        let distinct: HashSet<&str> = first.iter().map(String::as_str).collect();

        assert_eq!(first, second);
        assert_eq!(distinct.len(), first.len());
    }

    #[test]
    fn production_algorithms_agree_on_small_corpus() {
        let config = small_config();
        let mut results = Vec::new();
        for algorithm in Algorithm::ALL {
            let worker = WorkerConfig {
                k: 8,
                pattern_count: 16,
                algorithm,
                mode: SearchMode::All,
                output: PathBuf::new(),
            };
            let banks = build_banks(&config, &worker).unwrap();
            let order: Vec<usize> = (0..banks.len()).collect();
            results.push(search_banks(&banks, &order, SearchMode::All, 1));
        }

        assert!(results.windows(2).all(|pair| pair[0] == pair[1]));
    }
}
