#!/bin/bash

# Exit on error
set -e  

# Create patterns directory if it doesn't exist
mkdir -p ../patterns

# Function to generate reverse complement
reverse_complement() {
    local seq=$1
    echo "$seq" | tr 'ACGT' 'TGCA' | rev
}

# Function to generate k-mers from FASTQ
generate_kmers_fastq() {
    local input_file=$1
    local output_prefix=$2
    local k=$3
    
    echo "* Generating 100 ${k}-mers from FASTQ..."
    
    # Extract sequences from FASTQ and generate all k-mers
    echo "   - Generating k-mers..."
    awk -v k="$k" '
        NR % 4 == 2 {  # Only process sequence lines
            seq = $0
            if (length(seq) >= k) {  # Only process sequences that are long enough
                for (i = 1; i <= length(seq) - k + 1; i++) {
                    print substr(seq, i, k)
                }
            }
        }
    ' "$input_file" > "../patterns/${output_prefix}_${k}mers.tmp"
    
    echo "   - Selecting 100 random k-mers..."

    # Random source is seeded from repeatedly input string 'seed71'
    shuf --random-source=<(yes seed71) -n 100 "../patterns/${output_prefix}_${k}mers.tmp" | \
    awk '{print ">"NR"\n"$0}' > "../patterns/${output_prefix}_100x${k}mers.fasta"
    
    # Create grep-compatible version (without headers)
    grep -v "^>" "../patterns/${output_prefix}_100x${k}mers.fasta" > "../patterns/${output_prefix}_100x${k}mers.txt"
    
    # Create version with k-mers and their reverse complements
    echo "   - Generating k-mers with reverse complements..."
    while IFS= read -r kmer; do
        if [[ $kmer != ">"* ]]; then  # Skip FASTA headers
            echo "$kmer"
            reverse_complement "$kmer"
        fi
    done < "../patterns/${output_prefix}_100x${k}mers.fasta" > "../patterns/${output_prefix}_100x${k}mers_with_rc.txt"
    
    echo "   - Subsetting 1 ${k}-mer from FASTQ..."
    
    # Create single k-mer version using head
    head -2 "../patterns/${output_prefix}_100x${k}mers.fasta" > "../patterns/${output_prefix}_1x${k}mers.fasta"
    head -1 "../patterns/${output_prefix}_100x${k}mers.txt" > "../patterns/${output_prefix}_1x${k}mers.txt"
    
    # Create single k-mer version with reverse complement
    head -2 "../patterns/${output_prefix}_100x${k}mers_with_rc.txt" > "../patterns/${output_prefix}_1x${k}mers_with_rc.txt"
    
    # Clean up temporary file
    rm "../patterns/${output_prefix}_${k}mers.tmp"
}

# Generate k-mers from FASTQ
generate_kmers_fastq "../data/frag_1.fastq" "fastq" 31
generate_kmers_fastq "../data/frag_1.fastq" "fastq" 100

echo "K-mer generation complete!" 
