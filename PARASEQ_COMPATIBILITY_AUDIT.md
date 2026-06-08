# Paraseq Compatibility Audit

Phase: parser compatibility audit for migrating `extract` from `needletail` to `paraseq`.

Crate version audited: `paraseq 0.4.13`.

## Summary

`paraseq` is compatible with the core parser needs for `extract`:

- single-end FASTA parsing works;
- single-end FASTQ parsing works;
- paired-end FASTQ processing works through `fastx::Collection`;
- gz, bz2, and xz compressed FASTA inputs work through default `niffler` support;
- `Record::seq()` gives the correct newline-stripped sequence for matching;
- mismatched paired-end record-set sizes produce an error.

The main migration risk is output formatting:

- `needletail::SequenceRecord::write` preserves original record formatting more closely.
- `paraseq::Record::write_fasta` and `write_fastq` normalize output formatting.
- `paraseq::Record::seq_raw()` includes FASTA sequence line endings and should not be used for matching.

## API Findings

### Single-End Reader

`paraseq::fastx::Reader::from_path(path)` auto-detects FASTA/FASTQ and uses `niffler::send::from_path` when the default `niffler` feature is enabled.

Basic serial iteration:

```rust
let mut reader = paraseq::fastx::Reader::from_path(path)?;
let mut record_set = reader.new_record_set();

while record_set.fill(&mut reader)? {
    for record in record_set.iter() {
        let record = record?;
        let seq = record.seq();
    }
}
```

### Paired-End Reader

For two-file paired-end processing in `paraseq 0.4.13`, the clean API is `fastx::Collection`, not `reader1.process_parallel_paired(reader2, ...)`.

```rust
let collection = paraseq::fastx::Collection::from_paths(
    &[read_1, read_2],
    paraseq::fastx::CollectionType::Paired,
)?;

collection.process_parallel_paired(&mut processor, threads, threads_per_reader)?;
```

Migration consequence:

- Phase 6 should use `fastx::Collection` for paired-end migration.
- If we need lower-level manual chunking, we may need to inspect or wrap lower-level record-set APIs instead of only using built-in processors.

### Record Access

`paraseq::Record` exposes:

- `id() -> &[u8]`
- `seq() -> Cow<'_, [u8]>`
- `seq_raw() -> &[u8]`
- `qual() -> Option<&[u8]>`
- `write_fasta(...)`
- `write_fastq(...)`

Use `seq()` for pattern matching.

Do not use `seq_raw()` for matching because FASTA records can include line endings.

## Behavior Verified By Tests

Added integration tests in `tests/paraseq_compatibility.rs`.

Verified:

- `tests/fixtures/input/simple.fasta` parses as 3 FASTA records.
- `tests/fixtures/input/paired-1.fastq` parses as 2 FASTQ records with expected quality bytes.
- `tests/fixtures/input/fixed-width.faa` has:
  - `seq().len() == 280`
  - `seq_raw().len() == 284`
  - `seq_raw()` contains newline bytes
  - `seq()` does not contain newline bytes
- `tests/data/sample.fasta.gz`, `.bz2`, and `.xz` parse successfully and produce equivalent records.
- paired-end processing over `paired-1.fastq` and `paired-2.fastq` processes 2 pairs.
- mismatched paired-end lengths error with an incompatible record-set-size message.

## Compression

`paraseq` default features include `niffler` and `niffler/default`.

`niffler/default` includes:

- gzip;
- bgzip;
- bzip2;
- lzma/xz;
- zstd.

Migration consequence:

- Current compressed FASTA fixtures using `.gz`, `.bz2`, and `.xz` are supported.
- This is at least as broad as the current tested `needletail` compression coverage.

## Multiline FASTA

`paraseq::Record::seq()` strips newlines and returns `Cow<'_, [u8]>`.

Behavior:

- single-line FASTA can be borrowed;
- multiline FASTA requires newline removal and may allocate;
- `seq_raw()` avoids allocation but includes newlines and is not safe for matching.

Migration consequence:

- Use `seq()` for correctness.
- Accept possible allocation for multiline FASTA.
- For FASTQ, `seq()` is borrowed and should be efficient.

## Output Formatting Risk

`paraseq::Record::write_fasta` writes normalized FASTA:

```text
>id
seq
```

`paraseq::Record::write_fastq` writes normalized FASTQ:

```text
@id
seq
+
qual
```

Migration consequence:

- Existing output fixtures may change for fixed-width/multiline FASTA if we switch from `needletail` output writing to `paraseq` writing.
- This is probably acceptable if performance is the priority, but it should be treated as an intentional behavior decision.
- If exact original formatting is required, we need a different strategy, because `paraseq::Record` does not expose a `needletail`-style `all()` method in the public `Record` trait.

Recommendation:

- Normalize FASTA/FASTQ output in the `paraseq` path.
- Update docs/tests explicitly when Phase 5 changes output formatting.
- Preserve matching/logging semantics, not byte-for-byte multiline FASTA layout.

## Parallel Ordering Risk

The built-in `ParallelProcessor` and `PairedParallelProcessor` APIs expose records to workers but do not expose input-order record indices directly.

The internal implementation atomically claims batch positions, but the public processor callback only receives records or pairs.

Migration consequence:

- For deterministic ordered output/logging, do not rely on worker callback execution order.
- Options:
  - use a manual record-set pipeline with explicit chunk indices;
  - use `process_parallel_range` in deterministic ranges;
  - or build indexing into per-batch processing if the public API exposes enough batch order in practice.

Recommendation:

- Use serial `paraseq` first for Phase 5/6.
- For parallel `extract`, prefer a manual bounded record-set pipeline if deterministic output must remain byte-identical across runs.

## Paired-End Length Handling

`paraseq` detects incompatible paired record-set sizes and returns an error.

Migration consequence:

- Current MerKurio behavior for different paired-end file lengths can be preserved.
- Error wording will differ unless MerKurio maps the `paraseq` error into its current message.

Recommendation:

- Map paired-end length/record-set errors into MerKurio's existing user-facing error text.

## Phase 4/5 Recommendations

1. Add a parser-independent `FastxRecordView` and normalized writer.
2. Use `record.seq()` for matching.
3. Use `record.id()` for logging.
4. Use `record.qual()` only for FASTQ output.
5. Treat `seq_raw()` as audit/debug-only unless a future optimization needs raw FASTA access.
6. Start with serial `paraseq` single-end migration and update output-format expectations deliberately.
7. Use `fastx::Collection` for paired-end migration.
8. Defer built-in `process_parallel` usage until ordered aggregation is designed.

