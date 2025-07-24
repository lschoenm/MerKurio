# Minimal Example

Download the files `kmers.txt`, `sample.sam`, and `sample.fasta` from [the repository](https://github.com/lschoenm/MerKurio/example-minimal). 

After sucessfully installing MerKurio (follow the [official documentation](https://lschoenm.github.io/MerKurio/installation.html)) and having it added to your PATH, run the following commands to test MerKurio: 

```bash
# Run MerKurio extract
echo "Extracting sequences in sample.fasta containing the k-mer in kmers.txt:"
cmd="merkurio extract -f kmers.txt -i sample.fasta"
echo $cmd; eval $cmd

# Run MerKurio tag
echo "Tagging records in sample.sam with the k-mer in kmers.txt:"
cmd="merkurio tag -f kmers.txt -i sample.sam"
echo $cmd; eval $cmd

# Generate only a plain text log
echo "Generating a plain text log of matching statistics:"
cmd="merkurio tag -f kmers.txt -i sample.sam -S -l"
echo $cmd; eval $cmd

# Generate only a JSON log
echo "Generating a JSON log of matching statistics:"
cmd="merkurio tag -f kmers.txt -i sample.sam -S -j"
echo $cmd; eval $cmd
```

Explanation of the input files: 

- `kmers.txt`: contains a single _k_-mer sequence. 
- `sample.sam`: a synthetic example SAM file, containing three records. 
- `sample.fasta`: a synthetic FASTA file, containing three sequences. 

In this minimal example, no output files are generated. The records and matching statistics are directly written to the terminal. 