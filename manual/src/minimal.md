# Minimal Example

Download the files `kmers.txt`, `sample.sam`, and `sample.fasta` from [the repository](https://github.com/lschoenm/MerKurio/tree/master/example-minimal).

After sucessfully installing MerKurio (follow the [official documentation](https://lschoenm.github.io/MerKurio/installation.html)) and having it added to your PATH, run the following commands to test MerKurio:

```bash
# Run MerKurio extract
echo "Extracting sequences in sample.fasta containing the k-mer in kmers.txt:"
merkurio extract -f kmers.txt -i sample.fasta

# Run MerKurio tag
echo "Tagging records in sample.sam with the k-mer in kmers.txt:"
merkurio tag -f kmers.txt -i sample.sam

# Generate only a plain text log
echo "Generating a plain text log of matching statistics:"
merkurio tag -f kmers.txt -i sample.sam -S -l

# Generate only a JSON log
echo "Generating a JSON log of matching statistics:"
merkurio tag -f kmers.txt -i sample.sam -S -j
```

Explanation of the input files:

- `kmers.txt`: contains a single _k_-mer sequence.
- `sample.sam`: a synthetic example SAM file, containing three records.
- `sample.fasta`: a synthetic FASTA file, containing three sequences.

In this minimal example, no output files are generated. The records and matching statistics are directly written to the terminal.

In this minimal example, no output files are generated. The records and matching statistics are directly written to the terminal. For a description of the [plaint text log](./log.md) and [JSON log](./json.md), go to the official documentation.

For a detailed explanation of the available parameters of the [`tag`](./tag.md) and [`extract`](./extract.md) subcommands visit their sections in the documentation.
