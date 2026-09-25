#!/bin/bash

set -e  # Exit on error

# Create results directory if it doesn't exist
mkdir -p ../results

# Set number of warmup runs and benchmark runs
WARMUP=20
RUNS=100
# Maximum duration of each preflight run.
TIMEOUT=5m

# Replace with paths to the tested tools
MERKURIO="../../target/release/merkurio"
# Seqtool: https://github.com/markschl/seqtool/releases/tag/v0.4.0-beta.3
ST="../comparison_binaries/st-0.4.0-beta.3"
# Grep: https://www.gnu.org/software/grep/
GREP="../comparison_binaries/grep"
# fetch_reads: https://github.com/voichek/fetch_reads_with_kmers/releases/tag/V0_1_beta
FETCH_READS="../comparison_binaries/fetch_reads"
# Cookiecutter: https://github.com/ad3002/Cookiecutter/releases/tag/v1.0.0
CK="../comparison_binaries/cookiecutter-extract-1.0.0"
# SeqKit: https://github.com/shenwei356/seqkit/releases/tag/v2.10.0
SEQKIT="../comparison_binaries/seqkit"
# back_to_sequences: https://github.com/pierrepeterlongo/back_to_sequences
BACK_TO_SEQUENCES="back_to_sequences"

# Print system information and tool versions
echo "Machine running benchmarks:"
uname -a
echo "$(date)"
echo "Versions used:"
echo "* MerKurio: $($MERKURIO --version)"
echo "* seqtool: $($ST --version)"
echo "* fgrep: $($GREP --version)"
echo "* fetch_reads: no official release <https://github.com/voichek/fetch_reads_with_kmers>"
echo "* Cookiecutter 1.0.0"
echo "* seqkit: $($SEQKIT version)"
echo "* back_to_sequences: $($BACK_TO_SEQUENCES --version)"
echo ""
echo "Running benchmarks with $WARMUP warmup runs and $RUNS benchmark runs..."
echo "The 1,000,000 k-mer cases use 2 warmup runs and 10 benchmark runs."
echo "Preflight timeout: $TIMEOUT (forced termination after 5 additional seconds)"
echo "Using CPU 0 for all benchmarks (taskset -c 0) and highest priority (nice -20)"
echo ""

# Save versions to file
rm -f ../results/versions.txt
echo "Versions used ($(date)):" > ../results/versions.txt
echo "* MerKurio: $($MERKURIO --version)" >> ../results/versions.txt
echo "* seqtool: $($ST --version) <https://github.com/markschl/seqtool>" >> ../results/versions.txt
echo "* fgrep: $($GREP --version)" >> ../results/versions.txt
echo "* fetch_reads: no official release <https://github.com/voichek/fetch_reads_with_kmers>" >> ../results/versions.txt
echo "* Cookiecutter: 1.0.0" >> ../results/versions.txt
echo "* seqkit: $($SEQKIT version)" >> ../results/versions.txt
echo "* back_to_sequences: $($BACK_TO_SEQUENCES --version)" >> ../results/versions.txt

# Compare selected record IDs, ignoring order and extra header annotations.
record_ids() {
    awk -v format="$1" '
        (format == "fasta" && /^>/) || (format == "fastq" && NR % 4 == 1) {
            id = $1; sub(/^[@>]/, "", id); sub(/\r$/, "", id); print id
        }
    ' "$2" | LC_ALL=C sort
}

compare_outputs() {
    local format=$1 reference=$2
    shift 2
    local output
    for output in "$reference" "$@"; do
        if [[ ! -f "$output" ]]; then
            echo "Missing benchmark output: $output" >&2
            return 1
        fi
    done
    for output in "$@"; do
        if ! cmp -s <(record_ids "$format" "$reference") <(record_ids "$format" "$output"); then
            echo "Output mismatch: $output selects different record IDs than $reference" >&2
            return 1
        fi
    done
    echo "Output comparison passed: $reference"
}

preflight_benchmark() {
    local name=$1
    local command=$2
    local benchmark_command="nice -20 taskset -c 0 $command"
    BENCHMARK_ADDED=0

    echo "Checking $name..."
    if timeout --verbose --kill-after=5s "$TIMEOUT" bash -c "$benchmark_command"; then
        BENCHMARK_COMMANDS+=("$benchmark_command")
        BENCHMARK_ADDED=1
    elif [[ $? -eq 124 ]]; then
        echo "Skipping $name: exceeded the $TIMEOUT timeout"
    else
        echo "$name failed during the preflight run" >&2
        return 1
    fi
}

run_hyperfine() {
    local csv=$1
    local num_kmers=$2
    local warmup=$WARMUP
    local runs=$RUNS
    if [[ $num_kmers -eq 1000000 ]]; then
        warmup=2
        runs=10
    fi
    if [[ ${#BENCHMARK_COMMANDS[@]} -eq 0 ]]; then
        echo "No commands completed the preflight run"
        return
    fi
    rm -f "$csv"
    hyperfine --style color --warmup "$warmup" --runs "$runs" --export-csv "$csv" \
        "${BENCHMARK_COMMANDS[@]}"
}

# Function to run FASTA benchmarks
run_fasta_benchmarks() {
    local k=$1
    local num_kmers=$2
    local pattern_file="../patterns/fasta_${num_kmers}x${k}mers.fasta"
    local pattern_txt="../patterns/fasta_${num_kmers}x${k}mers.txt"
    local data_file="../data/genome-sl.fasta"
    local output_dir="../results/fasta"
    local st_ok=0 grep_ok=0 seqkit_ok=0 back_ok=0 merkurio_ok=0
    BENCHMARK_COMMANDS=()
    
    mkdir -p $output_dir
    
    echo -e "\n>>> Running benchmarks for ${num_kmers} x ${k} bp for FASTA"
    preflight_benchmark seqtool "$ST find -t 1 file:$pattern_file $data_file -o $output_dir/out-${num_kmers}x${k}mers-st.fasta -f"
    st_ok=$BENCHMARK_ADDED
    preflight_benchmark fgrep "$GREP -f $pattern_txt $data_file -B 1 --no-group-separator > $output_dir/out-${num_kmers}x${k}mers-fgrep.fasta"
    grep_ok=$BENCHMARK_ADDED
    preflight_benchmark seqkit "$SEQKIT grep -j 1 -P -s -f $pattern_txt $data_file > $output_dir/out-${num_kmers}x${k}mers-seqkit.fasta"
    seqkit_ok=$BENCHMARK_ADDED
    preflight_benchmark back_to_sequences "$BACK_TO_SEQUENCES --in-kmers $pattern_file --in-sequences $data_file --out-sequences $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fasta --out-kmers $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.kmers.fasta -k $k --stranded -t 1"
    back_ok=$BENCHMARK_ADDED
    preflight_benchmark MerKurio "$MERKURIO extract -i $data_file -f $pattern_file > $output_dir/out-${num_kmers}x${k}mers-merkurio.fasta"
    merkurio_ok=$BENCHMARK_ADDED

    run_hyperfine "$output_dir/${num_kmers}x${k}mers-results.csv" "$num_kmers"

    if (( merkurio_ok )); then
        local outputs=()
        (( st_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-st.fasta")
        (( grep_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-fgrep.fasta")
        (( seqkit_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-seqkit.fasta")
        (( back_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fasta")
        compare_outputs fasta "$output_dir/out-${num_kmers}x${k}mers-merkurio.fasta" "${outputs[@]}"
    else
        echo "Output comparison skipped because MerKurio timed out"
    fi
}

# Function to run FASTQ benchmarks
run_fastq_benchmarks() {
    local k=$1
    local num_kmers=$2
    local pattern_file="../patterns/fastq_${num_kmers}x${k}mers.fasta"
    local pattern_txt="../patterns/fastq_${num_kmers}x${k}mers.txt"
    local data_file="../data/frag_1.fastq"
    local data_file2="../data/frag_2.fastq"
    local output_dir="../results/fastq"
    local st_ok=0 grep_ok=0 ck_ok=0 seqkit_ok=0 back_ok=0 merkurio_ok=0
    BENCHMARK_COMMANDS=()
    
    mkdir -p $output_dir
    
    # Use `--seqtype other` to strictly match N characters
    echo -e "\n>>> Running benchmarks for ${num_kmers} x ${k} bp for FASTQ"
    preflight_benchmark seqtool "$ST find -t 1 file:$pattern_file $data_file -o $output_dir/out-${num_kmers}x${k}mers-st.fastq --seqtype other -f"
    st_ok=$BENCHMARK_ADDED
    preflight_benchmark fgrep "$GREP -f $pattern_txt $data_file -B 1 -A 2 --no-group-separator > $output_dir/out-${num_kmers}x${k}mers-fgrep.fastq"
    grep_ok=$BENCHMARK_ADDED
    preflight_benchmark Cookiecutter "$CK -i $data_file -f $pattern_txt -o $output_dir/out-${num_kmers}x${k}mers-ck"
    ck_ok=$BENCHMARK_ADDED
    preflight_benchmark seqkit "$SEQKIT grep -j 1 -P -s -f $pattern_txt $data_file > $output_dir/out-${num_kmers}x${k}mers-seqkit.fastq"
    seqkit_ok=$BENCHMARK_ADDED
    preflight_benchmark back_to_sequences "$BACK_TO_SEQUENCES --in-kmers $pattern_file --in-sequences $data_file --out-sequences $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fastq --out-kmers $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.kmers.fasta -k $k --stranded -t 1"
    back_ok=$BENCHMARK_ADDED
    preflight_benchmark MerKurio "$MERKURIO extract -i $data_file -f $pattern_file > $output_dir/out-${num_kmers}x${k}mers-merkurio.fastq"
    merkurio_ok=$BENCHMARK_ADDED

    run_hyperfine "$output_dir/${num_kmers}x${k}mers-results.csv" "$num_kmers"

    if (( merkurio_ok )); then
        local outputs=()
        (( st_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-st.fastq")
        (( grep_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-fgrep.fastq")
        (( ck_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-ck/frag_1.filtered.fastq")
        (( seqkit_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-seqkit.fastq")
        (( back_ok )) && outputs+=("$output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fastq")
        compare_outputs fastq "$output_dir/out-${num_kmers}x${k}mers-merkurio.fastq" "${outputs[@]}"
    else
        echo "Output comparison skipped because MerKurio timed out"
    fi
}

# Function to run paired-end FASTQ benchmarks (only for 31-mers) with reverse complements!
run_paired_end_benchmarks() {
    local num_kmers=$1
    local pattern_file="../patterns/fastq_${num_kmers}x31mers.fasta"
    local pattern_txt="../patterns/fastq_${num_kmers}x31mers.txt"
    local pattern_txt_rc="../patterns/fastq_${num_kmers}x31mers_with_rc.txt"
    local data_file="../data/frag_1.fastq"
    local data_file2="../data/frag_2.fastq"
    local output_dir="../results/fastq-paired"
    local fetch_ok=0 ck_ok=0 merkurio_ok=0
    BENCHMARK_COMMANDS=()
    
    mkdir -p $output_dir
    
    echo -e "\n>>> Running benchmarks for ${num_kmers} x 31 bp for paired-end FASTQ"
    preflight_benchmark fetch_reads "$FETCH_READS $data_file $data_file2 $pattern_file 31 $output_dir/out-${num_kmers}x31mers-fetch"
    fetch_ok=$BENCHMARK_ADDED
    preflight_benchmark Cookiecutter "$CK -1 $data_file -2 $data_file2 -f $pattern_txt_rc -o $output_dir/out-${num_kmers}x31mers-ck"
    ck_ok=$BENCHMARK_ADDED
    preflight_benchmark MerKurio "$MERKURIO extract -i $data_file -2 $data_file2 -f $pattern_file -o $output_dir/out-${num_kmers}x31mers -r"
    merkurio_ok=$BENCHMARK_ADDED

    run_hyperfine "$output_dir/${num_kmers}x31mers-results.csv" "$num_kmers"

    # Cookiecutter separates matching pairs and matching singleton reads.
    if (( ck_ok )); then
        cat "$output_dir/out-${num_kmers}x31mers-ck/frag_1.filtered.fastq" \
            "$output_dir/out-${num_kmers}x31mers-ck/frag_1.se.fastq" \
            > "$output_dir/out-${num_kmers}x31mers-ck-1.fastq"
        cat "$output_dir/out-${num_kmers}x31mers-ck/frag_2.filtered.fastq" \
            "$output_dir/out-${num_kmers}x31mers-ck/frag_2.se.fastq" \
            > "$output_dir/out-${num_kmers}x31mers-ck-2.fastq"
    fi

    if (( merkurio_ok )); then
        local outputs_1=() outputs_2=()
        (( fetch_ok )) && outputs_1+=("$output_dir/out-${num_kmers}x31mers-fetch_R1.fastq")
        (( fetch_ok )) && outputs_2+=("$output_dir/out-${num_kmers}x31mers-fetch_R2.fastq")
        (( ck_ok )) && outputs_1+=("$output_dir/out-${num_kmers}x31mers-ck-1.fastq")
        (( ck_ok )) && outputs_2+=("$output_dir/out-${num_kmers}x31mers-ck-2.fastq")
        compare_outputs fastq "$output_dir/out-${num_kmers}x31mers_1.fastq" "${outputs_1[@]}"
        compare_outputs fastq "$output_dir/out-${num_kmers}x31mers_2.fastq" "${outputs_2[@]}"
    else
        echo "Output comparison skipped because MerKurio timed out"
    fi
}


# Run FASTA benchmarks
# run_fasta_benchmarks 31 1
# run_fasta_benchmarks 31 100
# run_fasta_benchmarks 100 1
# run_fasta_benchmarks 100 100

# Run FASTQ benchmarks
run_fastq_benchmarks 31 1
run_fastq_benchmarks 31 100
run_fastq_benchmarks 31 1000000

# Run paired-end FASTQ benchmarks (only for 31-mers)
run_paired_end_benchmarks 1
run_paired_end_benchmarks 100
run_paired_end_benchmarks 1000000
echo "Benchmarks complete!"
