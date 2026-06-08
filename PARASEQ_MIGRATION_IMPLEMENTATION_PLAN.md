# Paraseq Migration Implementation Plan

## Goals

- Replace `needletail` in the `extract` subcommand with `paraseq`.
- Use `paraseq` record-set processing as the basis for true parallel pattern matching.
- Preserve current `extract` behavior for:
  - single-end FASTA/FASTQ;
  - paired-end FASTQ;
  - compressed inputs where supported;
  - output file naming;
  - plain text logs;
  - JSON logs;
  - `--invert-match`;
  - `--suppress-output`;
  - algorithm selection between Aho-Corasick and BNDMq.
- Preserve the fast non-verbose path that stops after the first match.
- Avoid unnecessary record copies where `paraseq` allows record-set-local borrowing.

## Non-Goals

- Do not migrate `tag`; it operates on SAM/BAM and is unrelated to FASTX parsing.
- Do not change matching algorithms during the parser migration.
- Do not require identical match ordering between Aho-Corasick and BNDMq unless current tests already require it.
- Do not add parallel log writing.
- Do not make `paraseq` the only implementation until behavior and benchmarks are validated.

## Rationale

The current `needletail` implementation is a sequential parser loop. Its records borrow parser-owned buffers, so sending records to worker threads requires copying IDs, sequences, and output bytes into owned work units.

`paraseq` is built around `RecordSet`s and parallel FASTX processing. That makes it a better fit for this problem because record sets are already the batching unit needed for parallelism, especially for paired-end files where synchronized record-set processing matters.

## Migration Strategy

Use a staged migration:

1. Introduce a common `PatternMatcher` API while still using `needletail`.
2. Add a serial `paraseq` implementation for single-end extract.
3. Add a serial `paraseq` implementation for paired-end extract.
4. Add parallel `paraseq` processing.
5. Benchmark `needletail` serial, `paraseq` serial, and `paraseq` parallel.
6. Remove or gate the old `needletail` extract path only after parity is proven.

This reduces risk by separating parser behavior changes from parallel matching changes.

## Phase 1: Matcher Abstraction

Add a parser-independent matching API.

Core types:

```rust
pub struct MatchHit {
    pub pattern_index: usize,
    pub position: usize,
}

pub enum PatternMatcher {
    AhoCorasick {
        ac: aho_corasick::AhoCorasick,
    },
    Bndmq {
        matchers: Vec<BNDMq>,
    },
}
```

Required methods:

```rust
impl PatternMatcher {
    pub fn new(
        pattern_list: &[String],
        use_aho_corasick: bool,
        case_insensitive: bool,
        q_size: Option<usize>,
    ) -> anyhow::Result<Self>;

    pub fn find_any(&self, seq: &[u8]) -> bool;

    pub fn for_each_match<F>(&self, seq: &[u8], on_match: F)
    where
        F: FnMut(MatchHit);

    pub fn algorithm_name(&self) -> &'static str;
}
```

Performance rules:

- `find_any` must remain allocation-free.
- `find_any` must preserve early exit for both algorithms.
- `for_each_match` must stream hits through a callback.
- Matching returns `pattern_index`, not cloned pattern strings.

Acceptance criteria:

- Existing `extract` tests pass with `needletail`.
- No performance regression in non-logging mode.

## Phase 2: Add `paraseq` Dependency

Update `Cargo.toml`.

Initial dependency:

```toml
paraseq = { version = "0.4", features = ["niffler"] }
```

Confirm feature names against the currently resolved crate before committing. If transparent decompression is available through a different feature or constructor, use that instead.

Keep `needletail` during migration:

```toml
needletail = { version = "0.7.1", features = ["compression"] }
```

Only remove `needletail` after all extract paths have moved and tests/benchmarks are acceptable.

## Phase 3: Parser Compatibility Audit

Before implementing, verify `paraseq` behavior for each current requirement.

Audit items:

- FASTQ single-end parsing.
- FASTA single-end parsing.
- Multiline FASTA behavior and whether `record.seq()` allocates.
- Gzip input support.
- Bzip2/xz input support if current `needletail` compression supports them.
- Record ID semantics compared with `needletail::SequenceRecord::id()`.
- Sequence semantics compared with `needletail::SequenceRecord::seq()`.
- Quality string access.
- Ability to reconstruct output records exactly enough for existing fixtures.
- Paired-end unequal-length error behavior.

Audit status:

- Completed in `PARASEQ_COMPATIBILITY_AUDIT.md`.
- Added guardrail tests in `tests/paraseq_compatibility.rs`.

Confirmed findings:

- `record.seq()` is the correct matching input.
- `record.seq_raw()` includes FASTA line endings and must not be used for matching.
- gz, bz2, and xz compressed FASTA inputs work through default `niffler` support.
- Two-file paired-end processing works through `fastx::Collection`.
- Paired-end unequal lengths produce a `paraseq` error that MerKurio should map to its current user-facing error text.

Expected issue:

- Exact output formatting may differ if `paraseq` reconstructs records instead of exposing original full record bytes. This must be tested against fixtures.

Decision point:

- If exact formatting preservation is required and `paraseq` cannot expose enough raw record data, either:
  - reconstruct output using MerKurio-owned formatting and update documented output behavior;
  - keep `needletail` for FASTA exact-format mode;
  - or support `paraseq` only for FASTQ parallel mode initially.

## Phase 4: Output Serialization Layer

Create a parser-independent output writer helper.

Status:

- Implemented in `src/fastx_output.rs`.
- The helper writes normalized FASTA/FASTQ records and validates FASTQ quality presence/length.
- Current `extract` still uses `needletail::SequenceRecord::write`; the normalized writer is intentionally not wired in until the `paraseq` path is implemented.

Suggested type:

```rust
enum FastxFormat {
    Fasta,
    Fastq,
}

struct FastxRecordView<'a> {
    id: &'a [u8],
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    format: FastxFormat,
}
```

Suggested helper:

```rust
fn write_fastx_record(
    writer: &mut dyn std::io::Write,
    record: FastxRecordView<'_>,
) -> anyhow::Result<()>;
```

Serialization policy:

- FASTQ output should be normalized to:

```text
@id
seq
+
qual
```

- FASTA output should be normalized to:

```text
>id
seq
```

Risk:

- Current `needletail` output preserves original line endings and multiline FASTA formatting. A normalized writer may change fixture output.

Recommendation:

- Use normalized output if performance is the priority and document this as acceptable.
- `paraseq::Record` does not expose a public `needletail`-style full-record `all()` method, so exact multiline FASTA layout preservation should not be assumed.

## Phase 5: Serial Single-End `paraseq` Extract

Implement a serial single-end path first.

Status:

- Implemented for single-end `extract`.
- Single-end input now uses `paraseq::fastx::Reader`.
- Single-end output now uses the normalized writer from `src/fastx_output.rs`.
- Paired-end input remains on `needletail` until Phase 6.
- The fixed-width FASTA output fixture was updated because multiline FASTA output is now normalized to one sequence line.

Flow:

```rust
let reader = paraseq::fastx::Reader::from_path(&args.in_fastx)?;
let mut record_set = reader.new_record_set();

while record_set.fill(&mut reader)? {
    for record in record_set.iter() {
        let record = record?;
        let seq = record.seq();

        if logging_active {
            matcher.for_each_match(seq, |hit| { ... });
        } else {
            matched = matcher.find_any(seq);
        }

        if matched != args.invert_match {
            write record unless suppressed;
        }
    }
}
```

Implementation notes:

- Avoid converting `seq` to `Vec<u8>` unless required by multiline FASTA normalization.
- Resolve record ID as borrowed bytes where possible.
- Update counters exactly as current code does.
- Keep writer and loggers on the main thread.

Acceptance criteria:

- All existing single-end `extract` tests pass or all intentional output-format differences are documented.
- `--suppress-output` works.
- Plain log and JSON log match current semantics.

## Phase 6: Serial Paired-End `paraseq` Extract

Implement paired-end processing using `fastx::Collection`.

Status:

- Implemented for paired-end `extract` with explicit synchronized `paraseq::fastx::Reader` record sets.
- The built-in `fastx::Collection`/processor API remains deferred because current output/log aggregation is serial and order-sensitive.
- Paired-end output now uses the normalized writer from `src/fastx_output.rs`.
- Added a direct `extract_records` test for mismatched paired-end file lengths.

Flow:

```rust
let mut reader_1 = paraseq::fastx::Reader::from_path(&args.in_fastx)?;
let mut reader_2 = paraseq::fastx::Reader::from_path(args.in_fastq_2.as_ref().unwrap())?;

let mut record_set_1 = reader_1.new_record_set();
let mut record_set_2 = reader_2.new_record_set();

// Fill both record sets, verify both advanced together, then zip records in order.
```

Rules:

- For non-logging mode:

```rust
matched = matcher.find_any(seq_1) || matcher.find_any(seq_2);
```

- For logging mode:

```rust
matcher.for_each_match(seq_1, |hit| { file_1 hit accounting });
matcher.for_each_match(seq_2, |hit| { file_2 hit accounting });
matched = file_1_matched || file_2_matched;
```

Acceptance criteria:

- Existing paired-end fixture passes.
- If file lengths differ, map the `paraseq` incompatible record-set-size error to MerKurio's current paired-end length error.
- Output files still use `_1` and `_2` suffixes.

## Phase 7: Parallel Processor Design

Use `paraseq`'s parallel processor traits for the parallel path.

Status:

- Added `src/extract_processing.rs` with owned worker-result types and pure record-processing helpers.
- Added `ExtractProcessor`, `SingleResult`, `PairedResult`, `RecordHit`, `RecordStats`, `OutputRecord`, and `RecordInput`.
- The helper preserves separate non-verbose first-hit and verbose all-hit paths.
- The helper does not write output/logs; it returns owned results suitable for Phase 8 ordered aggregation.
- Current `extract` command still uses the serial implementation directly until aggregation is introduced.

Shared processor state:

```rust
struct ExtractProcessor {
    matcher: Arc<PatternMatcher>,
    logging_active: bool,
    suppress_output: bool,
    invert_match: bool,
}
```

Per-thread processor state:

- A cloned processor should hold shared `Arc` data plus local buffers.
- Local buffers avoid lock contention.
- No worker should write to final output files or final log files.

Single-end worker result:

```rust
struct SingleResult {
    record_index: u64,
    matched: bool,
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
    hits: Vec<RecordHit>,
    stats: RecordStats,
}
```

Paired-end worker result:

```rust
struct PairedResult {
    pair_index: u64,
    matched: bool,
    read_1: OutputRecord,
    read_2: OutputRecord,
    hits: Vec<RecordHit>,
    stats_1: RecordStats,
    stats_2: RecordStats,
}
```

Important:

- Only include `id`, `seq`, and `qual` in results when output is not suppressed and the record should be extracted, or when logging needs `id`.
- In non-logging `--suppress-output` mode, result data can be very small: index, matched flag, and counters.
- In non-logging output mode, workers may need to copy output fields because the writer runs after worker callback returns.

## Phase 8: Ordered Aggregation

Parallel workers will complete out of order. Add an ordered aggregation layer.

Status:

- Added `OrderedResultBuffer<T>` in `src/extract_processing.rs`.
- Added `IndexedResult` implementations for single-end and paired-end worker results.
- Added `ExtractSummary` for merging record counts, base counts, hit counts, distinct-hit counts, extracted-record counts, and per-pattern counts.
- Added tests for out-of-order result draining and summary merging.
- Current `extract` command still uses direct serial aggregation; this layer is ready for the parallel pipeline.

Ordering mechanism:

- Assign monotonically increasing record or pair indices.
- Send results to an aggregator.
- Store out-of-order results in `BTreeMap<u64, Result>`.
- Drain while `next_expected_index` is present.

Aggregation responsibilities:

- Write records in input order.
- Write paired-end records in input order to separate output files.
- Emit plain logs in deterministic order.
- Emit JSON logs in deterministic order.
- Merge counters.
- Merge per-pattern hit counts.

Potential issue:

- `paraseq::ParallelProcessor::process_record` examples do not show index assignment. If the trait does not expose record-set or record indices, add indices during aggregation by processing record sets manually or using atomics carefully.

Indexing rule:

- Do not use `AtomicU64::fetch_add` inside workers unless it reflects input order. Worker execution order is not input order.
- Prefer indices derived from record-set order.

## Phase 9: Bounded Pipeline Alternative

If `paraseq`'s built-in parallel processor cannot preserve ordered output cleanly, use `paraseq` record sets manually with a bounded pipeline.

Status:

- Added `PipelineConfig` and `run_bounded_ordered_pipeline` in `src/extract_processing.rs`.
- The pipeline uses bounded work/result queues, producer and worker threads, and `OrderedResultBuffer` for deterministic draining.
- Added tests proving out-of-order worker completion is consumed in input order.
- Current `extract` command is not wired to this pipeline yet; Phase 10 will add thread control and the parallel path.

Pipeline:

```text
reader fills RecordSet chunks
  -> bounded work queue
  -> matcher workers process complete RecordSets
  -> bounded result queue
  -> ordered aggregator writes output/logs
```

Advantages:

- Explicit input-order chunk indices.
- Bounded memory.
- Clear writer ownership.
- Better control over copied output data.

Disadvantages:

- More implementation work than using `process_parallel`.
- Requires understanding `paraseq` record-set ownership APIs.

Recommendation:

- Try built-in `process_parallel` first for `--suppress-output` and logging-disabled output.
- Use manual record-set pipeline if ordered output/log preservation is awkward.

## Phase 10: CLI

Add thread control to `extract`.

Suggested option:

```text
-t, --threads <N>
```

Semantics:

- `--threads 1` uses serial `paraseq`.
- `--threads > 1` uses parallel `paraseq`.
- Default should remain `1` until benchmarks prove parallel mode is consistently beneficial.

Optional hidden option:

```text
--record-set-size <N>
```

Use this for benchmarking and tuning only.

## Phase 11: Tests

Parity tests:

- Existing single-end FASTA fixture.
- Existing single-end FASTA inverted fixture.
- Existing fixed-width FASTA fixture.
- Existing paired-end FASTQ fixture.
- Existing log-only fixture.
- Existing JSON fixtures.

Parallel determinism tests:

- `--threads 1` equals `--threads 2`.
- `--threads 1` equals `--threads 4`.
- Repeat `--threads 4` multiple times and verify byte-identical output/log/JSON.

Behavior tests:

- No matches.
- All records match.
- Dense matches with verbose logging.
- Sparse matches with verbose logging.
- `--suppress-output`.
- `--invert-match`.
- Paired-end mismatched file lengths.
- Aho-Corasick mode.
- BNDMq mode.
- Case-insensitive Aho-Corasick mode.

Compression tests:

- Gzip input.
- Bzip2 input if supported.
- Xz input if supported.

## Phase 12: Benchmarks

Benchmark these implementations:

- Current `needletail` serial.
- `paraseq` serial.
- `paraseq` parallel with 2, 4, 8, and available threads.

Benchmark dimensions:

- Single-end FASTQ.
- Paired-end FASTQ.
- FASTA.
- Compressed FASTQ.
- 1 pattern.
- 100 patterns.
- Long patterns.
- Sparse hits.
- Dense hits.
- Logging disabled.
- Plain log enabled.
- JSON log enabled.
- `--suppress-output`.

Expected outcomes:

- `paraseq` serial should be close to current `needletail` serial.
- `paraseq` parallel should help most with many patterns and output suppression.
- Logging-heavy workloads may remain writer-bound.
- Compressed workloads may remain decompression-bound.

## Phase 13: Rollout

Initial rollout:

- Keep `needletail` code behind a temporary internal fallback if migration risk remains.
- Add `--threads`; default to `1`.
- Use `paraseq` for `extract` in serial mode once parity is acceptable.
- Enable parallel mode only when `--threads > 1`.

After validation:

- Remove `needletail` from `extract`.
- Remove temporary fallback code.
- Consider removing `needletail` dependency entirely if no other code uses it.
- Update manual and benchmark documentation.

## Key Risks

- `paraseq` may not preserve exact output formatting currently produced by `needletail::SequenceRecord::write`.
- Parallel processor APIs may not expose input-order indices directly.
- Multiline FASTA may allocate during `seq()` extraction.
- Compression support may differ from `needletail`.
- Dense verbose logging may dominate runtime and hide parallel matching gains.
- `paraseq` is less established than `needletail`, so parser edge cases need tests.

## Risk Mitigations

- Keep migration staged and benchmarked.
- Preserve `find_any` as a distinct fast path.
- Keep writers/loggers single-threaded.
- Use ordered aggregation.
- Keep `needletail` fallback until parity is proven.
- Add fixture tests for exact output and log behavior.
- Benchmark before changing defaults.
