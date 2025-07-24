# Example Usage

An example usage of the `extract` subcommand to extract records from a FASTA file based on a list of _k_-mers and their reverse complements is shown below. Records from `input.fasta` matching any of the _k_-mers in `query_kmers.txt` are written to the file `output.fasta`. Logging information is written to the terminal (stdout):

```bash
merkurio extract -i input.fasta -o output.fasta -f query_kmers.txt -r -l
```

Another example where paired-end reads are extracted if they contain the sequence ACGT or TGCA; the extracted records are printed to the terminal and logging information is written to two files; as plain text to `log.txt`, and in JSON format to `log.json`:

```bash
merkurio extract -1 in_R1.fastq -2 in_R2.fastq -o output -s ACGT TGCA -l log.txt -j log.json
```

The following example only writes the matching statistics to a JSON file `log.json`; FASTA output is suppressed:

```bash
merkurio extract -i input.fasta -f query_kmers.txt -S -j log.json
```