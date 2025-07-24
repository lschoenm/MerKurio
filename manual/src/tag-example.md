# Example Usage

An example usage of the `tag` subcommand to tag a BAM file with the _k_-mers in the file `query_kmers.fasta` and output a SAM file is shown below:

```bash
merkurio tag -i input.bam -o output.sam -f query_kmers.fasta
```

Another example where the _k_-mers are provided on the command line and the search is also performed for their reverse complements. The tag is set to "MK". BAM file processing is done with 4 threads:

```bash
merkurio tag -i input.bam -o output.bam -s ACGT TGCA -r -p 4 -t MK
```

This example shows how to filter the output to only include records that contain at least one of the _k_-mers. Output is written to the terminal (stdout):

```bash
merkurio tag -i input.bam -s AAATCAAGA -m
```

The next example only writes the matching statistics to a JSON file `log.json`; SAM output is suppressed:

```bash
merkurio tag -i input.bam -s AAATCAAGA -S -j log.json
```