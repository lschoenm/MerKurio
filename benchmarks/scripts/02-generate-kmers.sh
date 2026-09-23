#!/bin/bash

# Exit on error
set -e  

# Create patterns directory if it doesn't exist
mkdir -p ../patterns

# Function to generate k-mers from FASTQ
generate_kmers_fastq() {
    local input_file=$1
    local output_prefix=$2
    local k=$3
    local count=${4:-100}
    
    echo "* Generating ${count} ${k}-mers from FASTQ..."
    
    # Extract sequences from FASTQ and generate all k-mers
    echo "   - Generating k-mers..."
    awk -v k="$k" '
        NR % 4 == 2 {  # Only process sequence lines
            seq = $0
            if (length(seq) >= k) {  # Only process sequences that are long enough
                for (i = 1; i <= length(seq) - k + 1; i++) {
                    kmer = substr(seq, i, k)
                    if (kmer !~ /[^ACGT]/) print kmer
                }
            }
        }
    ' "$input_file" > "../patterns/${output_prefix}_${k}mers.tmp"
    
    echo "   - Selecting ${count} random k-mers..."

    # Keep the million-query case at one million distinct patterns.
    if (( count == 1000000 )); then
        LC_ALL=C sort -u "../patterns/${output_prefix}_${k}mers.tmp" > "../patterns/${output_prefix}_${k}mers.unique.tmp"
        mv "../patterns/${output_prefix}_${k}mers.unique.tmp" "../patterns/${output_prefix}_${k}mers.tmp"
        if (( $(wc -l < "../patterns/${output_prefix}_${k}mers.tmp") < count )); then
            echo "Not enough distinct ${k}-mers for ${count} queries." >&2
            return 1
        fi
    fi

    # Random source is seeded from repeatedly input string 'seed71'
    shuf --random-source=<(yes seed71) -n "$count" "../patterns/${output_prefix}_${k}mers.tmp" | \
    awk '{print ">"NR"\n"$0}' > "../patterns/${output_prefix}_${count}x${k}mers.fasta"
    
    # Create grep-compatible version (without headers)
    grep -v "^>" "../patterns/${output_prefix}_${count}x${k}mers.fasta" > "../patterns/${output_prefix}_${count}x${k}mers.txt"
    
    # Create version with k-mers and their reverse complements
    echo "   - Generating k-mers with reverse complements..."
    awk '
        BEGIN { complement["A"]="T"; complement["T"]="A"; complement["C"]="G"; complement["G"]="C" }
        {
            print
            rc=""
            for (i=length($0); i>0; i--) {
                base=substr($0,i,1)
                rc=rc ((base in complement) ? complement[base] : base)
            }
            print rc
        }
    ' "../patterns/${output_prefix}_${count}x${k}mers.txt" > "../patterns/${output_prefix}_${count}x${k}mers_with_rc.txt"

    if (( count == 100 )); then
        echo "   - Subsetting 1 ${k}-mer from FASTQ..."
    
        # Create single k-mer version using head
        head -2 "../patterns/${output_prefix}_${count}x${k}mers.fasta" > "../patterns/${output_prefix}_1x${k}mers.fasta"
        head -1 "../patterns/${output_prefix}_${count}x${k}mers.txt" > "../patterns/${output_prefix}_1x${k}mers.txt"
    
        # Create single k-mer version with reverse complement
        head -2 "../patterns/${output_prefix}_${count}x${k}mers_with_rc.txt" > "../patterns/${output_prefix}_1x${k}mers_with_rc.txt"
    
    fi

    # Clean up temporary file
    rm "../patterns/${output_prefix}_${k}mers.tmp"
}

# Generate k-mers from FASTQ
generate_kmers_fastq "../data/frag_1.fastq" "fastq" 31
generate_kmers_fastq "../data/frag_1.fastq" "fastq" 100
generate_kmers_fastq "../data/frag_1.fastq" "fastq" 31 1000000

echo "K-mer generation complete!" 
