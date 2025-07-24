# Feature Overview

MerKurio provides two complementary subcommands:

- 🔍 **Extract**: Search FASTA/FASTQ data for _k_-mers and write records with matching _k_-mers to the terminal or a new file. 
  - Supports paired-end reads (a hit in one read extracts the whole pair).
- 📑 **Tag**: Annotate BAM/SAM alignments with _k_-mer tags and filter them based on matching _k_-mers.
  - Adds a two‑letter tag (default `km`) with comma-separated matching k‑mers (follows the [SAM format specification](https://samtools.github.io/hts-specs/SAMtags.pdf)). 
  - Optionally keeps only reads containing at least one _k_‑mer.
  - Multithreaded processing when working with BAM files.

Both commands share additional features: 

- Records [detailed matching statistics](./log.md) (positions of _k_-mer occurence, summary statistics, metadata). 
  - Human readable output in plain text. 
  - Structured JSON logs for easy machine parsing. 
- Reads compressed input files (`.gz`, `.bz2`, `.xz`). 
- Can seach for reverse complements or only canonical forms of _k_-mers. 
- Case-insensitive search or conversion to lower-/uppercase. 
- Inverse matching to keep only those records without matches. 
- Query _k_-mers can be provided as command line arguments or in a file (FASTA or plain text). 
- File types are inferred automatically. 
- Record output can be suppressed to only record statistics. 
