set -ex

# Run MerKurio extract
merkurio extract -f kmers.txt -i sample.fasta

echo ""

# Run MerKurio tag
merkurio tag -f kmers.txt -i sample.sam

echo ""

# Generate a plain text log
merkurio tag -f kmers.txt -i sample.sam -S -l

echo ""

# Generate a JSON log
merkurio tag -f kmers.txt -i sample.sam -S -j
