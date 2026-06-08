# Paraseq Native RecordSet Pipeline Plan

## Goal

Reduce parallel `extract` overhead by using `paraseq`'s native `RecordSet` buffers as pipeline work items instead of copying every input record into `OwnedRecordInput`.

The current parallel path is truly parallel, but it still pays a high memory cost before matching:

- The producer copies every record ID.
- The producer copies every sequence.
- The producer copies every FASTQ quality string.
- Most of those copies are wasted for non-matching records.
- In `--suppress-output`, all sequence/quality output copies are wasted.

The target behavior is:

- Match directly against borrowed records inside a `paraseq::fastx::RecordSet`.
- Copy only output records that are actually extracted.
- Copy record IDs only when logging needs them.
- Preserve deterministic output and log order.
- Keep the implementation readable enough to maintain.

## Why Not Use `process_parallel` Directly

`paraseq::fastx::Reader::process_parallel` and `process_parallel_paired` already process borrowed records without our per-record input copies. They are attractive for pure reductions, but they are not a clean fit for `extract` output.

Limitations for MerKurio:

- `process_parallel` does not expose stable batch start indices to the processor.
- Processor clones are thread-local, and `on_batch_complete` runs in worker completion order, not input order.
- Ordered FASTA/FASTQ output would need an additional ordering layer anyway.
- Ordered plain logs and JSON logs would need the same ordering layer.
- Writing from `on_batch_complete` through a mutex, as shown in `paraseq` examples, can reorder output batches.

Recommendation:

- Do not use `process_parallel` for normal `extract` output/log paths.
- Consider `process_parallel` only for a future stats-only mode where output/log ordering is irrelevant.
- Use manual `RecordSet` work chunks in our existing ordered bounded pipeline.

## Design

### Work Items

Replace owned record work chunks:

```rust
OwnedRecordInput
SingleWorkChunk { records: Vec<OwnedRecordInput> }
OwnedPairedInput
PairedWorkChunk { pairs: Vec<OwnedPairedInput> }
```

with native record-set chunks:

```rust
SingleRecordSetWorkChunk {
    start_index: u64,
    record_set: paraseq::fastx::RecordSet,
    format: FastxFormat,
}

PairedRecordSetWorkChunk {
    start_index: u64,
    record_set_1: paraseq::fastx::RecordSet,
    record_set_2: paraseq::fastx::RecordSet,
    format_1: FastxFormat,
    format_2: FastxFormat,
}
```

The producer owns the readers and fills record sets. Workers own the filled record sets while matching. The result consumer never sees borrowed data; it only receives owned output/log records already produced by `ExtractProcessor`.

### Processing

Workers iterate borrowed records from each `RecordSet`. For readability, a first implementation can use `record_set.iter()`:

```rust
for (offset, record) in chunk.record_set.iter().enumerate() {
    let record = record?;
    let seq = record.seq();
    processor.process_single_record(
        chunk.start_index + offset as u64,
        RecordInput {
            id: record.id(),
            seq: &seq,
            qual: record.qual(),
            format: chunk.format,
        },
    );
}
```

For the performance path, avoid `fastx::RecordSet::iter()` in the inner loop because it returns a boxed iterator. Dispatch once per chunk and then iterate over the concrete FASTA or FASTQ record set:

```rust
match &chunk.record_set {
    fastx::RecordSet::Fasta(records) => {
        for (offset, record) in records.iter().enumerate() {
            process_record(offset, record?, FastxFormat::Fasta)?;
        }
    }
    fastx::RecordSet::Fastq(records) => {
        for (offset, record) in records.iter().enumerate() {
            process_record(offset, record?, FastxFormat::Fastq)?;
        }
    }
}
```

This keeps dynamic dispatch out of per-record iteration while still localizing the format split to one helper.

This keeps all matching borrowed from `paraseq`'s record-set buffer. `ExtractProcessor` already converts records to owned `OutputRecord` only when output is needed. It already converts hit IDs to owned bytes only when logging creates `RecordHit`s.

For paired-end input, compare record-set lengths before processing and then zip records:

```rust
let len_1 = record_set_len(&chunk.record_set_1);
let len_2 = record_set_len(&chunk.record_set_2);
if len_1 != len_2 {
    bail!("The two input files have a different number of records...");
}
```

Then process each pair with `pair_index = chunk.start_index + offset`.

## Implementation Steps

1. Add record-set work chunk types.

Create `SingleRecordSetWorkChunk` and `PairedRecordSetWorkChunk` in `src/cmd_extract.rs`. Add a small `record_set_len(&fastx::RecordSet) -> usize` helper by matching `RecordSet::Fasta` and `RecordSet::Fastq` and calling `n_records()`.

2. Add record-set processing functions.

Add `process_single_record_set_chunk` and `process_paired_record_set_chunk`. These should call the existing `ExtractProcessor` with borrowed `RecordInput`s and return `Vec<SingleResult>` / `Vec<PairedResult>`.

Use one small helper for single-record processing and one for paired-record processing so the concrete FASTA/FASTQ match arms do not duplicate extraction logic.

3. Replace parallel producers.

In the `worker_threads > 1` single-end branch, replace manual `OwnedRecordInput` chunk filling with:

- `reader.new_record_set_with_size(EXTRACT_PARALLEL_CHUNK_SIZE)`
- `record_set.fill(&mut reader)`
- send `SingleRecordSetWorkChunk { start_index, record_set, format }`
- increment `start_index` by `record_set_len`
- allocate the next record set after the previous one has been moved into the work queue

In the paired-end branch, fill two record sets together and send one `PairedRecordSetWorkChunk`. Keep the existing mismatched-file-length error semantics.

4. Keep the ordered result pipeline.

Keep `run_bounded_ordered_pipeline_with_producer` and the chunked result channel. It already reduces result-channel traffic and preserves deterministic output. The work item changes should be transparent to the consumer.

5. Remove owned input structs.

Delete `OwnedRecordInput`, `SingleWorkChunk`, `OwnedPairedInput`, `PairedWorkChunk`, `process_single_work_chunk`, and `process_paired_work_chunk` after the record-set path is working.

6. Add optional record-set pooling only if needed.

Initial implementation can allocate one `RecordSet` per chunk. That is already much cheaper than copying every record. If benchmarks still show allocation pressure, add a bounded pool:

- Producer receives empty record sets from a pool channel.
- Producer fills them and sends them to workers.
- Workers return processed record sets to the pool after processing.
- `RecordSet::fill()` clears and reuses the existing internal buffers, so returned sets can be refilled without per-chunk buffer reallocation.
- Pool size should be about `worker_count + work_queue_bound`.

Do not add pooling before measuring unless allocation shows up clearly; it adds complexity and error-path handling.

## Performance Considerations

Expected wins:

- Removes all producer-side ID/sequence/quality copies for non-matching records.
- Removes quality copying entirely in `--suppress-output`.
- In sparse-hit output mode, copies only the small fraction of records that are extracted.
- Keeps matching close to `paraseq`'s native batch buffers.
- Reduces allocator pressure substantially compared with `OwnedRecordInput` per record.
- Avoids boxed iterator overhead in the hot loop if the concrete per-format iteration helper is used.

Costs that remain:

- FASTQ parsing is still bounded by reader throughput.
- Output/log writing is still single-consumer and ordered.
- Extracted records still need owned copies before the record-set buffer can be dropped.
- JSON/plain logging still copies record IDs for each hit.
- `record.seq()` may allocate for multiline FASTA because normalized sequence may need to be materialized.
- Creating one fresh record set per chunk still allocates metadata and may allocate input buffers; add pooling only if benchmarks or profiling show this is material.

Chunk size:

- Keep `EXTRACT_PARALLEL_CHUNK_SIZE = 8192` initially.
- For 150 bp paired FASTQ, one paired chunk is roughly 5 MB of input buffer data.
- With bounded queues, memory should remain predictable.
- Benchmark `2048`, `8192`, and `16384` if performance remains inconclusive.

Thread count:

- Benchmark explicit `--threads 2` and `--threads 4`.
- Avoid interpreting `--threads 0` as the main scaling result because it may use too many logical cores and increase memory contention.

Pattern workload:

- `250` x `31 bp` patterns may still be too cheap for parallelism if Aho-Corasick dominates efficiently.
- Also benchmark `1000`, `5000`, and `10000` patterns.
- Use both sparse-hit and dense-hit inputs because output/log costs change the bottleneck.

## Validation

Correctness tests:

- Existing serial-vs-parallel single-end fixture parity.
- Existing serial-vs-parallel paired-end fixture parity.
- Parallel determinism repeated runs.
- Paired-end mismatched file length error.
- Fixed-width FASTA sequence normalization.
- JSON and plain log pattern counts.
- `--suppress-output`.

Performance benchmarks:

- `RECORDS=1000000 PATTERNS=250 THREADS="1 2 4" RUNS=3`
- `RECORDS=1000000 PATTERNS=1000 THREADS="1 2 4" RUNS=3`
- `RECORDS=1000000 PATTERNS=5000 THREADS="1 2 4" RUNS=3`
- Compare `fast-output`, `json-output`, and `json-suppress` separately.

Success criteria:

- Parallel mode should reduce wall time in `fast-output` for sufficiently high pattern counts.
- Parallel mode should not regress correctness or output order.
- Memory use should be lower than the owned-input pipeline for large paired FASTQ.
- If wall time still does not improve, profile reader throughput and matcher CPU time before adding more pipeline complexity.

## Risk Assessment

Main risks:

- Borrowed record data must not escape worker processing.
- Ordered output must remain byte-identical to serial output.
- Paired record-set length mismatches must keep the current error behavior.
- Record-set pooling can deadlock if the consumer exits early and workers try to return buffers.

Mitigation:

- Keep result data owned, as it is today.
- Keep ordered consumer unchanged.
- Implement without pooling first.
- Add pooling only after benchmark evidence justifies it.
