#!/bin/bash

set -e  # Exit on error

# Create data directory if it doesn't exist
mkdir -p ../data

# Download FASTA file
echo "Downloading FASTA file..."
curl -L -o ../data/genome.fasta https://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta

# Download FASTQ files
echo "Downloading FASTQ files..."
curl -L -o ../data/frag_1.fastq.gz https://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_1.fastq.gz
curl -L -o ../data/frag_2.fastq.gz https://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/frag_2.fastq.gz

# Decompress FASTQ files
echo "Decompressing FASTQ files..."
gunzip -f ../data/frag_1.fastq.gz
gunzip -f ../data/frag_2.fastq.gz

# Remove compressed FASTQ files
echo "Removing compressed FASTQ files..."
rm -f ../data/frag_1.fastq.gz
rm -f ../data/frag_2.fastq.gz

# Generate single-line FASTA file
echo "Generating single-line FASTA file..."
awk '/^>/{if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' ../data/genome.fasta > ../data/genome-sl.fasta

echo "Download complete!" 