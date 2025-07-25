# Quick Example

After sucessfully installing MerKurio (follow the [official documentation](https://lschoenm.github.io/MerKurio/installation.html)) and having it added to your PATH, download the files `kmers.txt`, `sample.sam`, and `sample.fasta`. Then run the script `example.sh` or the following commands to test MerKurio:

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

In this minimal example, no output files are generated. The records and matching statistics are directly written to the terminal. For a description of the [plaint text log](https://lschoenm.github.io/MerKurio/log.html) and [JSON log](https://lschoenm.github.io/MerKurio/json.html), go to the official documentation.

For a detailed explanation of the available parameters of the [`tag`](https://lschoenm.github.io/MerKurio/tag.html) and [`extract`](https://lschoenm.github.io/MerKurio/extract.html) subcommands visit the documentation.
