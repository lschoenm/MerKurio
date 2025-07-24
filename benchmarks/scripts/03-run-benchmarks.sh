#!/bin/bash

set -e  # Exit on error

# Create results directory if it doesn't exist
mkdir -p ../results

# Set number of warmup runs and benchmark runs
WARMUP=20
RUNS=100

# Replace with paths to the tested tools
MERKURIO="../../target/release/merkurio"
# Seqtool: https://github.com/markschl/seqtool/releases/tag/v0.4.0-beta.3
ST="../st-0.4.0-beta.3"
# Grep: https://www.gnu.org/software/grep/
GREP="../grep"
# fetch_reads: https://github.com/voichek/fetch_reads_with_kmers/releases/tag/V0_1_beta
FETCH_READS="../fetch_reads"
# Cookiecutter: https://github.com/ad3002/Cookiecutter/releases/tag/v1.0.0
CK="../cookiecutter-extract-1.0.0"
# SeqKit: https://github.com/shenwei356/seqkit/releases/tag/v2.10.0
SEQKIT="../seqkit"
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
echo "Using CPU 0 for all benchmarks (taskset -c 0) and highest priority (nice -20)"
echo ""

# Save versions to file
echo "Versions used ($(date)):" > ../results/versions.txt
echo "* MerKurio: $($MERKURIO --version)" >> ../results/versions.txt
echo "* seqtool: $($ST --version) <https://github.com/markschl/seqtool>" >> ../results/versions.txt
echo "* fgrep: $($GREP --version)" >> ../results/versions.txt
echo "* fetch_reads: no official release <https://github.com/voichek/fetch_reads_with_kmers>" >> ../results/versions.txt
echo "* Cookiecutter: 1.0.0" >> ../results/versions.txt
echo "* seqkit: $($SEQKIT version)" >> ../results/versions.txt
echo "* back_to_sequences: $($BACK_TO_SEQUENCES --version)" >> ../results/versions.txt

# Function to run FASTA benchmarks
run_fasta_benchmarks() {
    local k=$1
    local num_kmers=$2
    local pattern_file="../patterns/fasta_${num_kmers}x${k}mers.fasta"
    local pattern_txt="../patterns/fasta_${num_kmers}x${k}mers.txt"
    local data_file="../data/genome-sl.fasta"
    local output_dir="../results/fasta"
    
    mkdir -p $output_dir
    
    echo -e "\n>>> Running benchmarks for ${num_kmers} x ${k} bp for FASTA"
    hyperfine --style color --warmup $WARMUP --runs $RUNS --export-csv $output_dir/${num_kmers}x${k}mers-results.csv \
        "nice -20 taskset -c 0 $ST find -t 1 file:$pattern_file $data_file -o $output_dir/out-${num_kmers}x${k}mers-st.fasta -f" \
        "nice -20 taskset -c 0 $GREP -f $pattern_txt $data_file -B 1 --no-group-separator > $output_dir/out-${num_kmers}x${k}mers-fgrep.fasta" \
        "nice -20 taskset -c 0 $SEQKIT grep -j 1 -P -s -f $pattern_txt $data_file > $output_dir/out-${num_kmers}x${k}mers-seqkit.fasta" \
        "nice -20 taskset -c 0 $BACK_TO_SEQUENCES --in-kmers $pattern_file --in-sequences $data_file --out-sequences $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fasta --out-kmers $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.kmers.fasta -k $k --stranded -t 1" \
        "nice -20 taskset -c 0 $MERKURIO extract -i $data_file -f $pattern_file > $output_dir/out-${num_kmers}x${k}mers-merkurio.fasta"
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
    
    mkdir -p $output_dir
    
    # Use `--seqtype other` to strictly match N characters
    echo -e "\n>>> Running benchmarks for ${num_kmers} x ${k} bp for FASTQ"
    hyperfine --style color --warmup $WARMUP --runs $RUNS --export-csv $output_dir/${num_kmers}x${k}mers-results.csv \
        "nice -20 taskset -c 0 $ST find -t 1 file:$pattern_file $data_file -o $output_dir/out-${num_kmers}x${k}mers-st.fastq --seqtype other -f" \
        "nice -20 taskset -c 0 $GREP -f $pattern_txt $data_file -B 1 -A 2 --no-group-separator > $output_dir/out-${num_kmers}x${k}mers-fgrep.fastq" \
        "nice -20 taskset -c 0 $CK -i $data_file -f $pattern_txt -o $output_dir/out-${num_kmers}x${k}mers-ck" \
        "nice -20 taskset -c 0 $SEQKIT grep -j 1 -P -s -f $pattern_txt $data_file > $output_dir/out-${num_kmers}x${k}mers-seqkit.fastq" \
        "nice -20 taskset -c 0 $BACK_TO_SEQUENCES --in-kmers $pattern_file --in-sequences $data_file --out-sequences $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.fastq --out-kmers $output_dir/out-${num_kmers}x${k}mers-back_to_sequences.kmers.fasta -k $k --stranded -t 1" \
        "nice -20 taskset -c 0 $MERKURIO extract -i $data_file -f $pattern_file > $output_dir/out-${num_kmers}x${k}mers-merkurio.fastq"
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
    
    mkdir -p $output_dir
    
    echo -e "\n>>> Running benchmarks for ${num_kmers} x 31 bp for paired-end FASTQ"
    hyperfine --style color --warmup $WARMUP --runs $RUNS --export-csv $output_dir/${num_kmers}x31mers-results.csv \
        "nice -20 taskset -c 0 $FETCH_READS $data_file $data_file2 $pattern_file 31 $output_dir/out-${num_kmers}x31mers-fetch" \
        "nice -20 taskset -c 0 $CK -1 $data_file -2 $data_file2 -f $pattern_txt_rc -o $output_dir/out-${num_kmers}x31mers-ck" \
        "nice -20 taskset -c 0 $MERKURIO extract -i $data_file -2 $data_file2 -f $pattern_file -o $output_dir/out-${num_kmers}x31mers -r"
}


# Run FASTA benchmarks
# run_fasta_benchmarks 31 1
# run_fasta_benchmarks 31 100
# run_fasta_benchmarks 100 1
# run_fasta_benchmarks 100 100

# Run FASTQ benchmarks
run_fastq_benchmarks 31 1
run_fastq_benchmarks 31 100

# Run paired-end FASTQ benchmarks (only for 31-mers)
run_paired_end_benchmarks 1
run_paired_end_benchmarks 100
echo "Benchmarks complete!" 