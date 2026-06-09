#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd "$BENCH_DIR/.." && pwd)"
RESULT_DIR="$BENCH_DIR/results/realistic-extract"
DATA_DIR="$BENCH_DIR/data/realistic-extract"
MERKURIO="${MERKURIO:-$REPO_DIR/target/release/merkurio}"
MERKURIO_SINGLE="${MERKURIO_SINGLE:-$BENCH_DIR/merkurio-single}"
RUNS="${RUNS:-5}"
WARMUP="${WARMUP:-1}"
RECORDS="${RECORDS:-1000000}"
READ_LEN="${READ_LEN:-150}"
PATTERNS="${PATTERNS:-250}"
HIT_EVERY="${HIT_EVERY:-50}"
THREADS="${THREADS:-1 2 4}"
CHUNK_SIZES="${CHUNK_SIZES:-}"

mkdir -p "$RESULT_DIR" "$DATA_DIR"

if [[ ! -x "$MERKURIO" ]]; then
    echo "Missing MerKurio release binary at $MERKURIO" >&2
    echo "Build it first with: cargo build --release" >&2
    exit 1
fi

DATA_TAG="${RECORDS}x${READ_LEN}_${PATTERNS}p_hit${HIT_EVERY}"
FASTQ_1="$DATA_DIR/realistic_${DATA_TAG}_1.fastq"
FASTQ_2="$DATA_DIR/realistic_${DATA_TAG}_2.fastq"
PATTERN_FILE="$DATA_DIR/patterns_${PATTERNS}x31.txt"
METADATA_FILE="$DATA_DIR/metadata_${DATA_TAG}.txt"

generate_inputs() {
    if [[ -s "$FASTQ_1" && -s "$FASTQ_2" && -s "$PATTERN_FILE" && -s "$METADATA_FILE" ]]; then
        return
    fi

    python3 - "$FASTQ_1" "$FASTQ_2" "$PATTERN_FILE" "$METADATA_FILE" "$RECORDS" "$READ_LEN" "$PATTERNS" "$HIT_EVERY" <<'PY'
import random
import sys

fastq_1, fastq_2, pattern_file, metadata_file = sys.argv[1:5]
record_count = int(sys.argv[5])
read_len = int(sys.argv[6])
pattern_count = int(sys.argv[7])
hit_every = int(sys.argv[8])

if read_len < 31:
    raise SystemExit("READ_LEN must be at least 31")
if pattern_count < 1:
    raise SystemExit("PATTERNS must be at least 1")
if hit_every < 1:
    raise SystemExit("HIT_EVERY must be at least 1")

rng = random.Random(20260608)
bases = "ACGT"

def rand_seq(length):
    return "".join(rng.choice(bases) for _ in range(length))

patterns = []
seen = set()
while len(patterns) < pattern_count:
    pattern = rand_seq(31)
    if pattern not in seen:
        seen.add(pattern)
        patterns.append(pattern)

with open(pattern_file, "w", encoding="ascii") as handle:
    for pattern in patterns:
        handle.write(f"{pattern}\n")

hit_records_r1 = 0
hit_records_r2 = 0
total_bases = record_count * read_len * 2
qual = "I" * read_len

with open(fastq_1, "w", encoding="ascii") as r1, open(fastq_2, "w", encoding="ascii") as r2:
    for index in range(record_count):
        seq_1 = rand_seq(read_len)
        seq_2 = rand_seq(read_len)

        if index % hit_every == 0:
            pattern = patterns[index % pattern_count]
            pos = (index * 17) % (read_len - len(pattern) + 1)
            seq_1 = seq_1[:pos] + pattern + seq_1[pos + len(pattern):]
            hit_records_r1 += 1

        if index % (hit_every * 3) == 0:
            pattern = patterns[(index * 7) % pattern_count]
            pos = (index * 11) % (read_len - len(pattern) + 1)
            seq_2 = seq_2[:pos] + pattern + seq_2[pos + len(pattern):]
            hit_records_r2 += 1

        r1.write(f"@read_{index}/1\n{seq_1}\n+\n{qual}\n")
        r2.write(f"@read_{index}/2\n{seq_2}\n+\n{qual}\n")

with open(metadata_file, "w", encoding="ascii") as handle:
    handle.write(f"records={record_count}\n")
    handle.write(f"read_len={read_len}\n")
    handle.write(f"patterns={pattern_count}\n")
    handle.write("pattern_len=31\n")
    handle.write(f"hit_every={hit_every}\n")
    handle.write(f"expected_r1_hit_records={hit_records_r1}\n")
    handle.write(f"expected_r2_hit_records={hit_records_r2}\n")
    handle.write(f"total_bases={total_bases}\n")
PY
}

run_hyperfine_group() {
    local mode="$1"
    shift
    local output_csv="$RESULT_DIR/realistic-extract-${DATA_TAG}-${mode}.csv"

    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-csv "$output_csv" \
        "$@"
    echo "Wrote $output_csv"
}

run_shell_timer_group() {
    local mode="$1"
    shift
    local output="$RESULT_DIR/summary-${DATA_TAG}-${mode}.txt"
    : > "$output"

    for command in "$@"; do
        echo ">>> $command" | tee -a "$output"
        for run in $(seq 1 "$RUNS"); do
            /usr/bin/time -p bash -c "$command" 2>&1 | awk -v run="$run" '/^real / { print "run " run ": " $2 " s" }' | tee -a "$output"
        done
        echo "" | tee -a "$output"
    done
    echo "Wrote $output"
}

run_benchmark_group() {
    local mode="$1"
    shift

    if command -v hyperfine >/dev/null 2>&1; then
        run_hyperfine_group "$mode" "$@"
    else
        echo "hyperfine not found; using /usr/bin/time fallback for $mode."
        run_shell_timer_group "$mode" "$@"
    fi
}

main() {
    generate_inputs

    echo "Realistic extract parallel benchmarks"
    echo "Binary: $MERKURIO"
    if [[ -x "$MERKURIO_SINGLE" ]]; then
        echo "Single-threaded baseline: $MERKURIO_SINGLE"
    else
        echo "Single-threaded baseline: not found at $MERKURIO_SINGLE; skipping baseline"
    fi
    echo "Records: $RECORDS read pairs"
    echo "Read length: $READ_LEN"
    echo "Patterns: $PATTERNS x 31 bp"
    echo "Hit every: $HIT_EVERY records in R1, every $((HIT_EVERY * 3)) records in R2"
    echo "Threads: $THREADS"
    if [[ -n "$CHUNK_SIZES" ]]; then
        echo "Chunk sizes: $CHUNK_SIZES"
    else
        echo "Chunk sizes: default"
    fi
    echo "Runs: $RUNS, warmup: $WARMUP"
    echo "Results: $RESULT_DIR"

    local fast_output_commands=()
    if [[ -x "$MERKURIO_SINGLE" ]]; then
        fast_output_commands+=("$MERKURIO_SINGLE extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE -o $RESULT_DIR/fast-output-${DATA_TAG}-single")
    fi
    if [[ -n "$CHUNK_SIZES" ]]; then
        for chunk_size in $CHUNK_SIZES; do
            for thread_count in $THREADS; do
                fast_output_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count --chunk-size $chunk_size -o $RESULT_DIR/fast-output-${DATA_TAG}-t${thread_count}-c${chunk_size}")
            done
        done
    else
        for thread_count in $THREADS; do
            fast_output_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count -o $RESULT_DIR/fast-output-${DATA_TAG}-t${thread_count}")
        done
    fi

    local json_output_commands=()
    if [[ -x "$MERKURIO_SINGLE" ]]; then
        json_output_commands+=("$MERKURIO_SINGLE extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE -o $RESULT_DIR/json-output-${DATA_TAG}-single -j $RESULT_DIR/json-output-${DATA_TAG}-single.json")
    fi
    if [[ -n "$CHUNK_SIZES" ]]; then
        for chunk_size in $CHUNK_SIZES; do
            for thread_count in $THREADS; do
                json_output_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count --chunk-size $chunk_size -o $RESULT_DIR/json-output-${DATA_TAG}-t${thread_count}-c${chunk_size} -j $RESULT_DIR/json-output-${DATA_TAG}-t${thread_count}-c${chunk_size}.json")
            done
        done
    else
        for thread_count in $THREADS; do
            json_output_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count -o $RESULT_DIR/json-output-${DATA_TAG}-t${thread_count} -j $RESULT_DIR/json-output-${DATA_TAG}-t${thread_count}.json")
        done
    fi

    local json_suppress_commands=()
    if [[ -x "$MERKURIO_SINGLE" ]]; then
        json_suppress_commands+=("$MERKURIO_SINGLE extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --suppress-output -j $RESULT_DIR/json-suppress-${DATA_TAG}-single.json")
    fi
    if [[ -n "$CHUNK_SIZES" ]]; then
        for chunk_size in $CHUNK_SIZES; do
            for thread_count in $THREADS; do
                json_suppress_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count --chunk-size $chunk_size --suppress-output -j $RESULT_DIR/json-suppress-${DATA_TAG}-t${thread_count}-c${chunk_size}.json")
            done
        done
    else
        for thread_count in $THREADS; do
            json_suppress_commands+=("$MERKURIO extract -i $FASTQ_1 -2 $FASTQ_2 -f $PATTERN_FILE --threads $thread_count --suppress-output -j $RESULT_DIR/json-suppress-${DATA_TAG}-t${thread_count}.json")
        done
    fi

    run_benchmark_group "fast-output" "${fast_output_commands[@]}"
    run_benchmark_group "json-output" "${json_output_commands[@]}"
    run_benchmark_group "json-suppress" "${json_suppress_commands[@]}"
}

main "$@"
