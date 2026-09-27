#![allow(dead_code)]

use std::collections::HashSet;
use std::env;
use std::fs::{self, File};
use std::hint::black_box;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[path = "../../../src/pattern_matching.rs"]
mod pattern_matching;
#[path = "../../../src/pattern_preprocessing.rs"]
mod pattern_preprocessing;

use pattern_matching::BNDMq;

const DEFAULT_SEED: u64 = 0x4d65_724b_7572_696f;

#[derive(Clone, Copy)]
enum SearchMode {
    First,
    All,
}

impl SearchMode {
    fn name(self) -> &'static str {
        match self {
            Self::First => "first",
            Self::All => "all",
        }
    }
}

struct Config {
    k_min: usize,
    k_max: usize,
    max_q: usize,
    sequence_length: usize,
    sequences_per_pattern: usize,
    patterns_per_k: usize,
    runs: usize,
    target_run_ms: u64,
    seed: u64,
    output: PathBuf,
    metadata: PathBuf,
}

impl Default for Config {
    fn default() -> Self {
        let results = Path::new(env!("CARGO_MANIFEST_DIR")).join("results");
        Self {
            k_min: 4,
            k_max: usize::BITS.min(64) as usize,
            max_q: 12,
            sequence_length: 150,
            sequences_per_pattern: 512,
            patterns_per_k: 256,
            runs: 15,
            target_run_ms: 20,
            seed: DEFAULT_SEED,
            output: results.join("q_sweep.csv"),
            metadata: results.join("metadata.txt"),
        }
    }
}

impl Config {
    fn parse() -> Result<Self, String> {
        let mut config = Self::default();
        let mut args = env::args().skip(1);
        while let Some(flag) = args.next() {
            let value = |args: &mut std::iter::Skip<std::env::Args>,
                         flag: &str|
             -> Result<String, String> {
                args.next()
                    .ok_or_else(|| format!("Missing value for {flag}."))
            };
            match flag.as_str() {
                "--k-min" => config.k_min = parse_usize(&value(&mut args, &flag)?, &flag)?,
                "--k-max" => config.k_max = parse_usize(&value(&mut args, &flag)?, &flag)?,
                "--max-q" => config.max_q = parse_usize(&value(&mut args, &flag)?, &flag)?,
                "--sequence-length" => {
                    config.sequence_length = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--sequences" => {
                    config.sequences_per_pattern = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--patterns-per-k" => {
                    config.patterns_per_k = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--runs" | "--samples" => {
                    config.runs = parse_usize(&value(&mut args, &flag)?, &flag)?
                }
                "--target-ms" => {
                    config.target_run_ms = parse_u64(&value(&mut args, &flag)?, &flag)?
                }
                "--seed" => config.seed = parse_u64(&value(&mut args, &flag)?, &flag)?,
                "--output" => config.output = PathBuf::from(value(&mut args, &flag)?),
                "--metadata" => config.metadata = PathBuf::from(value(&mut args, &flag)?),
                "-h" | "--help" => {
                    print_help();
                    std::process::exit(0);
                }
                _ => return Err(format!("Unknown option: {flag}. Use --help for usage.")),
            }
        }
        config.validate()?;
        Ok(config)
    }

    fn validate(&self) -> Result<(), String> {
        let word_bits = usize::BITS as usize;
        if self.k_min < 1 || self.k_min > self.k_max {
            return Err("Require 1 <= k-min <= k-max.".to_string());
        }
        if self.k_max > word_bits {
            return Err(format!(
                "BNDMq supports at most {word_bits}-base patterns on this architecture."
            ));
        }
        if self.k_max > self.sequence_length {
            return Err("k-max cannot exceed sequence-length.".to_string());
        }
        if self.max_q == 0
            || self.sequences_per_pattern == 0
            || self.patterns_per_k == 0
            || self.runs == 0
            || self.target_run_ms == 0
        {
            return Err(
                "max-q, sequences, patterns-per-k, runs, and target-ms must be positive."
                    .to_string(),
            );
        }
        if self.sequences_per_pattern % 2 != 0 {
            return Err(
                "sequences must be even so the mixed corpus has exactly 50% hits.".to_string(),
            );
        }
        let maximum_patterns = 4_usize.checked_pow(self.k_min as u32).unwrap_or(usize::MAX);
        if self.patterns_per_k > maximum_patterns {
            return Err(format!(
                "patterns-per-k exceeds the {maximum_patterns} distinct DNA patterns possible at k-min={}.",
                self.k_min
            ));
        }
        Ok(())
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

fn print_help() {
    println!(
        "\
MerKurio BNDMq q-value tuning

Options:
  --k-min N             Smallest pattern length [default: 4]
  --k-max N             Largest pattern length [default: machine word size]
  --max-q N             Largest q to test [default: 12]
  --sequence-length N   Bases per sequence [default: 150]
  --sequences N         Sequences per pattern [default: 512]
  --patterns-per-k N    Distinct patterns for each k [default: 256]
  --runs N              Measurement rounds [default: 15]
  --target-ms N         Approximate minimum duration per timed run [default: 20]
  --seed N              Deterministic unsigned integer seed
  --output PATH         Raw CSV output path
  --metadata PATH       Run metadata output path
  -h, --help            Show this help"
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

fn contains_pattern(sequence: &[u8], pattern: &[u8]) -> bool {
    sequence
        .windows(pattern.len())
        .any(|window| window == pattern)
}

fn naive_matches(sequence: &[u8], pattern: &[u8]) -> Vec<usize> {
    sequence
        .windows(pattern.len())
        .enumerate()
        .filter_map(|(position, window)| (window == pattern).then_some(position))
        .collect()
}

fn no_hit_sequence(rng: &mut SplitMix64, length: usize, pattern: &[u8]) -> Vec<u8> {
    loop {
        let sequence = random_dna(rng, length);
        if !contains_pattern(&sequence, pattern) {
            return sequence;
        }
    }
}

struct PatternCase {
    pattern: Vec<u8>,
    mixed_corpus: Vec<Vec<u8>>,
    expected_matches: Vec<Vec<usize>>,
}

fn build_pattern_bank(
    seed: u64,
    k: usize,
    sequence_length: usize,
    sequences_per_pattern: usize,
    pattern_count: usize,
) -> Vec<PatternCase> {
    let mut rng = SplitMix64::new(seed ^ (k as u64).wrapping_mul(0xd6e8_feb8_6659_fd93));
    let mut patterns = HashSet::with_capacity(pattern_count);
    let mut cases = Vec::with_capacity(pattern_count);

    while cases.len() < pattern_count {
        let pattern = random_dna(&mut rng, k);
        if !patterns.insert(pattern.clone()) {
            continue;
        }
        let mut mixed_corpus: Vec<Vec<u8>> = (0..sequences_per_pattern)
            .map(|_| no_hit_sequence(&mut rng, sequence_length, &pattern))
            .collect();
        for (index, sequence) in mixed_corpus.iter_mut().enumerate() {
            if index % 2 == 0 {
                let position = rng.range(sequence_length - k + 1);
                sequence[position..position + k].copy_from_slice(&pattern);
            }
        }
        let expected_matches = mixed_corpus
            .iter()
            .map(|sequence| naive_matches(sequence, &pattern))
            .collect();
        cases.push(PatternCase {
            pattern,
            mixed_corpus,
            expected_matches,
        });
    }
    cases
}

fn build_matchers(pattern_bank: &[PatternCase], q: usize) -> Result<Vec<BNDMq>, String> {
    pattern_bank
        .iter()
        .map(|case| BNDMq::new(&case.pattern, q).map_err(|error| error.to_string()))
        .collect()
}

fn validate_matchers(
    matchers: &[BNDMq],
    pattern_bank: &[PatternCase],
    k: usize,
    q: usize,
) -> Result<(), String> {
    for (pattern_index, (matcher, case)) in matchers.iter().zip(pattern_bank).enumerate() {
        for (sequence_index, (sequence, expected)) in case
            .mixed_corpus
            .iter()
            .zip(&case.expected_matches)
            .enumerate()
        {
            let observed: Vec<usize> = matcher.find_iter(sequence).collect();
            if &observed != expected {
                return Err(format!(
                    "Incorrect all-match result at k={k}, q={q}, pattern={pattern_index}, sequence={sequence_index}."
                ));
            }
            if matcher.find_match(sequence) != !expected.is_empty() {
                return Err(format!(
                    "Incorrect first-match result at k={k}, q={q}, pattern={pattern_index}, sequence={sequence_index}."
                ));
            }
        }
    }
    Ok(())
}

fn search_pattern_bank(
    matchers: &[BNDMq],
    pattern_bank: &[PatternCase],
    pattern_order: &[usize],
    mode: SearchMode,
    iterations: usize,
) -> u64 {
    let mut checksum = 0_u64;
    for _ in 0..iterations {
        for &pattern_index in pattern_order {
            let matcher = &matchers[pattern_index];
            for sequence in &pattern_bank[pattern_index].mixed_corpus {
                let sequence = black_box(sequence.as_slice());
                match mode {
                    SearchMode::First => {
                        checksum = checksum.wrapping_add(u64::from(matcher.find_match(sequence)));
                    }
                    SearchMode::All => {
                        for position in matcher.find_iter(sequence) {
                            checksum = checksum.wrapping_add(position as u64 + 1);
                        }
                    }
                }
            }
        }
    }
    black_box(checksum)
}

fn timed_search(
    matchers: &[BNDMq],
    pattern_bank: &[PatternCase],
    pattern_order: &[usize],
    mode: SearchMode,
    iterations: usize,
) -> (u128, u64) {
    let start = Instant::now();
    let checksum = search_pattern_bank(matchers, pattern_bank, pattern_order, mode, iterations);
    (start.elapsed().as_nanos(), checksum)
}

fn calibrated_iterations(
    matchers: &[BNDMq],
    pattern_bank: &[PatternCase],
    pattern_order: &[usize],
    mode: SearchMode,
    target_run_ms: u64,
) -> usize {
    let (elapsed_ns, _) = timed_search(matchers, pattern_bank, pattern_order, mode, 1);
    let target_ns = u128::from(target_run_ms) * 1_000_000;
    target_ns.div_ceil(elapsed_ns.max(1)).clamp(1, 10_000) as usize
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
        command_output("git", &["status", "--short"]).replace('\n', "; ")
    )?;
    writeln!(
        output,
        "rustc={}",
        command_output("rustc", &["--version", "--verbose"])
    )?;
    writeln!(output, "system={}", command_output("uname", &["-a"]))?;
    writeln!(output, "corpus=mixed_50_percent_hits")?;
    writeln!(output, "k_min={}", config.k_min)?;
    writeln!(output, "k_max={}", config.k_max)?;
    writeln!(output, "max_q={}", config.max_q)?;
    writeln!(output, "sequence_length={}", config.sequence_length)?;
    writeln!(
        output,
        "sequences_per_pattern={}",
        config.sequences_per_pattern
    )?;
    writeln!(output, "patterns_per_k={}", config.patterns_per_k)?;
    writeln!(output, "runs={}", config.runs)?;
    writeln!(output, "target_run_ms={}", config.target_run_ms)?;
    writeln!(output, "seed={}", config.seed)?;
    Ok(())
}

fn run(config: &Config) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(parent) = config.output.parent() {
        fs::create_dir_all(parent)?;
    }
    write_metadata(config)?;
    let mut output = BufWriter::new(File::create(&config.output)?);
    writeln!(
        output,
        "k,q,mode,run,iterations,patterns,sequences_per_pattern,sequence_length,bases,elapsed_ns,ns_per_base,checksum"
    )?;

    for k in config.k_min..=config.k_max {
        eprintln!("tuning k={k}/{}: generating pattern bank", config.k_max);
        let pattern_bank = build_pattern_bank(
            config.seed,
            k,
            config.sequence_length,
            config.sequences_per_pattern,
            config.patterns_per_k,
        );
        let natural_pattern_order: Vec<usize> = (0..pattern_bank.len()).collect();
        let q_values: Vec<usize> = (1..=config.max_q.min(k)).collect();
        let mut matcher_banks = Vec::with_capacity(q_values.len());
        for &q in &q_values {
            let matchers = build_matchers(&pattern_bank, q)?;
            validate_matchers(&matchers, &pattern_bank, k, q)?;
            matcher_banks.push(matchers);
        }

        for mode in [SearchMode::First, SearchMode::All] {
            eprintln!("tuning k={k}/{} mode={}", config.k_max, mode.name());
            let iterations: Vec<usize> = matcher_banks
                .iter()
                .map(|matchers| {
                    calibrated_iterations(
                        matchers,
                        &pattern_bank,
                        &natural_pattern_order,
                        mode,
                        config.target_run_ms,
                    )
                })
                .collect();
            for (q_index, matchers) in matcher_banks.iter().enumerate() {
                black_box(search_pattern_bank(
                    matchers,
                    &pattern_bank,
                    &natural_pattern_order,
                    mode,
                    iterations[q_index],
                ));
            }

            let mut round_rng = SplitMix64::new(
                config.seed
                    ^ (k as u64).wrapping_mul(0xa076_1d64_78bd_642f)
                    ^ match mode {
                        SearchMode::First => 0xe703_7ed1_a0b4_28db,
                        SearchMode::All => 0x8ebc_6af0_9c88_c6e3,
                    },
            );
            for run_index in 0..config.runs {
                let mut q_order: Vec<usize> = (0..q_values.len()).collect();
                let mut pattern_order = natural_pattern_order.clone();
                round_rng.shuffle(&mut q_order);
                round_rng.shuffle(&mut pattern_order);

                for q_index in q_order {
                    let iteration_count = iterations[q_index];
                    let (elapsed_ns, checksum) = timed_search(
                        &matcher_banks[q_index],
                        &pattern_bank,
                        &pattern_order,
                        mode,
                        iteration_count,
                    );
                    let bases = iteration_count
                        .saturating_mul(config.patterns_per_k)
                        .saturating_mul(config.sequences_per_pattern)
                        .saturating_mul(config.sequence_length);
                    let ns_per_base = elapsed_ns as f64 / bases as f64;
                    writeln!(
                        output,
                        "{k},{},{},{run_index},{iteration_count},{},{},{},{bases},{elapsed_ns},{ns_per_base:.9},{checksum}",
                        q_values[q_index],
                        mode.name(),
                        config.patterns_per_k,
                        config.sequences_per_pattern,
                        config.sequence_length,
                    )?;
                }
            }
        }
        output.flush()?;
    }
    Ok(())
}

fn main() {
    let config = Config::parse().unwrap_or_else(|error| {
        eprintln!("Error: {error}");
        std::process::exit(2);
    });
    if let Err(error) = run(&config) {
        eprintln!("Tuning failed: {error:#}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tuning_tests {
    use super::*;

    #[test]
    fn generated_pattern_bank_is_distinct_and_half_matching() {
        let cases = build_pattern_bank(DEFAULT_SEED, 4, 150, 32, 16);
        let distinct: HashSet<&[u8]> = cases.iter().map(|case| case.pattern.as_slice()).collect();

        assert_eq!(distinct.len(), 16);
        for case in cases {
            assert_eq!(
                case.expected_matches
                    .iter()
                    .filter(|matches| !matches.is_empty())
                    .count(),
                16
            );
        }
    }

    #[test]
    fn shuffle_is_deterministic_and_complete() {
        let mut first: Vec<usize> = (1..=12).collect();
        let mut second = first.clone();
        SplitMix64::new(DEFAULT_SEED).shuffle(&mut first);
        SplitMix64::new(DEFAULT_SEED).shuffle(&mut second);
        let mut sorted = first.clone();
        sorted.sort_unstable();

        assert_eq!(first, second);
        assert_eq!(sorted, (1..=12).collect::<Vec<_>>());
    }
}
