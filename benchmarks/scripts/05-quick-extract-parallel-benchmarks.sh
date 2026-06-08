#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd "$BENCH_DIR/.." && pwd)"
RESULT_DIR="$BENCH_DIR/results/quick-extract"
DATA_DIR="$BENCH_DIR/data/quick-extract"
MERKURIO="${MERKURIO:-$REPO_DIR/target/release/merkurio}"
RUNS="${RUNS:-5}"
WARMUP="${WARMUP:-1}"

mkdir -p "$RESULT_DIR" "$DATA_DIR"

if [[ ! -x "$MERKURIO" ]]; then
    echo "Missing MerKurio release binary at $MERKURIO" >&2
    echo "Build it first with: cargo build --release" >&2
    exit 1
fi

generate_inputs() {
    local fasta="$DATA_DIR/quick.fasta"
    local fastq_1="$DATA_DIR/quick_1.fastq"
    local fastq_2="$DATA_DIR/quick_2.fastq"
    local patterns="$DATA_DIR/patterns.txt"
    local many_patterns="$DATA_DIR/patterns-many.txt"

    if [[ -s "$fasta" && -s "$fastq_1" && -s "$fastq_2" && -s "$patterns" && -s "$many_patterns" ]]; then
        return
    fi

    python3 - "$fasta" "$fastq_1" "$fastq_2" "$patterns" "$many_patterns" <<'PY'
import sys

fasta, fastq_1, fastq_2, patterns, many_patterns = sys.argv[1:]
motif = "ACGTTGCA"
miss = "TTTTGGGGCCCCAAAA"
record_count = 20_000

with open(patterns, "w", encoding="ascii") as handle:
    handle.write(f"{motif}\n")

with open(many_patterns, "w", encoding="ascii") as handle:
    for i in range(100):
        handle.write(("ACGT" + f"{i:04d}" + "TGCA").replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T").replace("4", "A").replace("5", "C").replace("6", "G").replace("7", "T").replace("8", "A").replace("9", "C") + "\n")
    handle.write(f"{motif}\n")

with open(fasta, "w", encoding="ascii") as handle:
    for i in range(record_count):
        seq = (motif if i % 10 == 0 else miss) * 8
        handle.write(f">record_{i}\n{seq}\n")

with open(fastq_1, "w", encoding="ascii") as r1, open(fastq_2, "w", encoding="ascii") as r2:
    for i in range(record_count):
        seq_1 = (motif if i % 10 == 0 else miss) * 8
        seq_2 = (motif if i % 13 == 0 else miss) * 8
        r1.write(f"@record_{i}/1\n{seq_1}\n+\n{'I' * len(seq_1)}\n")
        r2.write(f"@record_{i}/2\n{seq_2}\n+\n{'I' * len(seq_2)}\n")
PY
}

run_hyperfine() {
    local name="$1"
    shift
    local output_csv="$RESULT_DIR/${name}.csv"

    hyperfine \
        --warmup "$WARMUP" \
        --runs "$RUNS" \
        --export-csv "$output_csv" \
        "$@"
}

run_shell_timer() {
    local output="$RESULT_DIR/summary.txt"
    : > "$output"

    for command in "$@"; do
        echo ">>> $command" | tee -a "$output"
        for run in $(seq 1 "$RUNS"); do
            /usr/bin/time -p bash -c "$command" 2>&1 | awk -v run="$run" '/^real / { print "run " run ": " $2 " s" }' | tee -a "$output"
        done
        echo "" | tee -a "$output"
    done
}

main() {
    generate_inputs

    local fasta="$DATA_DIR/quick.fasta"
    local fastq_1="$DATA_DIR/quick_1.fastq"
    local fastq_2="$DATA_DIR/quick_2.fastq"
    local patterns="$DATA_DIR/patterns.txt"
    local many_patterns="$DATA_DIR/patterns-many.txt"

    local single_output="$RESULT_DIR/single.fasta"
    local paired_output="$RESULT_DIR/paired"

    local commands=(
        "$MERKURIO extract -i $fasta -f $patterns --threads 1 --suppress-output -j $RESULT_DIR/single-t1.json"
        "$MERKURIO extract -i $fasta -f $patterns --threads 2 --suppress-output -j $RESULT_DIR/single-t2.json"
        "$MERKURIO extract -i $fasta -f $patterns --threads 4 --suppress-output -j $RESULT_DIR/single-t4.json"
        "$MERKURIO extract -i $fasta -f $patterns --threads 0 --suppress-output -j $RESULT_DIR/single-t0.json"
        "$MERKURIO extract -i $fasta -f $many_patterns --threads 1 -o $single_output"
        "$MERKURIO extract -i $fasta -f $many_patterns --threads 4 -o $single_output"
        "$MERKURIO extract -i $fastq_1 -2 $fastq_2 -f $patterns --threads 1 --suppress-output -j $RESULT_DIR/paired-t1.json"
        "$MERKURIO extract -i $fastq_1 -2 $fastq_2 -f $patterns --threads 4 --suppress-output -j $RESULT_DIR/paired-t4.json"
        "$MERKURIO extract -i $fastq_1 -2 $fastq_2 -f $patterns --threads 0 --suppress-output -j $RESULT_DIR/paired-t0.json"
        "$MERKURIO extract -i $fastq_1 -2 $fastq_2 -f $many_patterns --threads 4 -o $paired_output"
    )

    echo "Quick extract benchmarks"
    echo "Binary: $MERKURIO"
    echo "Runs: $RUNS, warmup: $WARMUP"
    echo "Results: $RESULT_DIR"

    if command -v hyperfine >/dev/null 2>&1; then
        run_hyperfine "quick-extract" "${commands[@]}"
        echo "Wrote $RESULT_DIR/quick-extract.csv"
    else
        echo "hyperfine not found; using /usr/bin/time fallback."
        run_shell_timer "${commands[@]}"
        echo "Wrote $RESULT_DIR/summary.txt"
    fi
}

main "$@"
