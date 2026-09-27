#!/bin/bash

set -e

# Run from the scripts directory, like 03-run-benchmarks.sh.
WARMUP=20
RUNS=100
TIMEOUT=5m  # Preflight only
MERKURIO="../../target/release/merkurio"
ST="../comparison_binaries/st-0.4.0-beta.3"
SEQKIT="../comparison_binaries/seqkit"
BACK_TO_SEQUENCES="back_to_sequences"

# Separate from results so script 04 does not include this additional check.
OUTPUT_DIR="../results-multithreaded"
mkdir -p "$OUTPUT_DIR"

# One logical CPU per physical core, using the topology reported by Linux.
mapfile -t CORES < <(lscpu -p=CPU,CORE,SOCKET | awk -F, '!/^#/ && !seen[$3 ":" $2]++ { print $1 }')
if (( ${#CORES[@]} < 6 )); then
    echo "This check needs at least six physical cores" >&2
    exit 1
fi

uname -a
date
"$MERKURIO" --version
"$ST" --version
"$SEQKIT" version
"$BACK_TO_SEQUENCES" --version

# Detect the output format: back_to_sequences 0.8.4 writes FASTA even for FASTQ input.
record_ids() {
    awk '
        NR == 1 { fasta = /^>/ }
        (fasta && /^>/) || (!fasta && NR % 4 == 1) {
            id = $1; sub(/^[@>]/, "", id); sub(/\r$/, "", id); print id
        }
    ' "$1" | LC_ALL=C sort
}

preflight_benchmark() {
    local name=$1 output=$2
    local command="nice -20 taskset -c $cpus $3"
    local status
    echo "Checking $name..."
    if timeout --verbose --kill-after=5s "$TIMEOUT" bash -c "$command"; then
        commands+=("$command")
        outputs+=("$output")
    else
        status=$?
        if [[ $status -eq 124 || $status -eq 137 ]]; then
            echo "Skipping $name: preflight timed out or was killed (exit $status)"
        else
            echo "$name failed during preflight (exit $status)" >&2
            return "$status"
        fi
    fi
}

run_benchmark() {
    local threads=$1 num_kmers=$2
    local mode=${3:-comparison}
    local cpus
    cpus=$(IFS=,; echo "${CORES[*]:0:threads}")
    local pattern="../patterns/fastq_${num_kmers}x31mers.fasta"
    local pattern_txt="../patterns/fastq_${num_kmers}x31mers.txt"
    local reads="../data/frag_1.fastq"
    local prefix="$OUTPUT_DIR/${num_kmers}x31mers-${threads}threads"
    local warmup=$WARMUP runs=$RUNS
    local commands=() outputs=()
    if [[ $mode == scaling || $num_kmers -eq 1000000 ]]; then
        warmup=2
        runs=10
    fi

    if [[ $mode == scaling ]]; then
        prefix="$OUTPUT_DIR/${num_kmers}x31mers-merkurio-scaling"
        echo -e "\n>>> MerKurio: $num_kmers x 31 bp, 1/2/4/6 threads"
    else
        echo -e "\n>>> $num_kmers x 31 bp, $threads threads, CPUs $cpus"
    fi
    echo "$warmup warmups, $runs measured runs; preflight timeout $TIMEOUT"
    # Remove the old CSV even if every tool is skipped this time.
    rm -f "$prefix-results.csv"
    if [[ $mode == scaling ]]; then
        for threads in 1 2 4 6; do
            cpus=$(IFS=,; echo "${CORES[*]:0:threads}")
            preflight_benchmark "MerKurio ($threads threads, CPUs $cpus)" "$prefix-${threads}threads.fastq" \
                "$MERKURIO extract -i $reads -f $pattern --threads $threads > $prefix-${threads}threads.fastq"
        done
    else
        preflight_benchmark MerKurio "$prefix-merkurio.fastq" \
            "$MERKURIO extract -i $reads -f $pattern --threads $threads > $prefix-merkurio.fastq"
        preflight_benchmark seqtool "$prefix-st.fastq" \
            "$ST find -t $threads file:$pattern $reads -o $prefix-st.fastq --seqtype other -f"
        preflight_benchmark seqkit "$prefix-seqkit.fastq" \
            "$SEQKIT grep -j $threads -P -s -f $pattern_txt $reads > $prefix-seqkit.fastq"
        preflight_benchmark back_to_sequences "$prefix-back_to_sequences.fastq" \
            "$BACK_TO_SEQUENCES --in-kmers $pattern --in-sequences $reads --out-sequences $prefix-back_to_sequences.fastq -k 31 --stranded -t $threads"
    fi

    if (( ${#commands[@]} == 0 )); then
        echo "No commands completed preflight"
        return
    fi
    hyperfine --style color --warmup "$warmup" --runs "$runs" \
        --export-csv "$prefix-results.csv" "${commands[@]}"

    # Compare successful tools against the first successful output.
    local output
    for output in "${outputs[@]}"; do
        if [[ ! -f "$output" ]]; then
            echo "Missing benchmark output: $output" >&2
            return 1
        fi
    done
    for output in "${outputs[@]:1}"; do
        if ! cmp -s <(record_ids "${outputs[0]}") <(record_ids "$output"); then
            echo "Output mismatch: $output differs from ${outputs[0]}" >&2
            return 1
        fi
    done
}

# Compare MerKurio thread counts before comparing programs.
for num_kmers in 1 100 1000000; do
    run_benchmark 1 "$num_kmers" scaling
done

for num_kmers in 1 100 1000000; do
    run_benchmark 4 "$num_kmers"
done
echo "Multithreaded check complete! Results: $OUTPUT_DIR"
