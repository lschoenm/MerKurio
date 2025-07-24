# The `tag` Subcommand

Running `merkurio tag` will tag aligned sequences in a [SAM/BAM file](https://doi.org/10.1093/bioinformatics/btp352) with contained _k_-mers. If a record contains one or more of the _k_-mers, it is annotated with a tag ("km" by default; must be exactly two characters long) and the respective _k_-mers, separated by commas. 

The _k_-mers can be provided as a list of strings on the command line or in a file (FASTA or plain text file). The output is written to a SAM or BAM file, depending on the file extension. If the chosen tag is already present in the records, the new values are added to existing ones. Optionally, keep only records which are matching at least one _k_-mer, or keep only records which are not matching any _k_-mer. Multithreading is supported for parsing BAM files. 

Note that searching with MerKurio **is case-sensitive** by default, but can be set to ignore the case. Alternatively, all input _k_-mers can be converted to lower- or uppercase. 

**Detailed match statistics** are written to the terminal (stdout) or to a file if a path is provided, showing which records got matched by which patterns, including the zero-based start of the position of the match. Matching statistics can also be saved in [**JSON**](https://developer.mozilla.org/en-US/docs/Learn/JavaScript/Objects/JSON) **format** for parsing by other programs. 

Further options include the ability to also search for the **reverse complements** of nucleotide query sequences (i.e., a sequence is reversed, and C/G and A/T are swapped), or searching for canonical _k_-mers only. Output of matching records can be suppressed if only the statistics are needed. MerKurio tries to select the most efficient algorithm for the given query sequences, but the user can override this choice by selecting a specific algorithm (not recommended).

For a documentation of the SAM format, see the [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf). \
For a documentation on the SAM/BAM tags, see the [SAM optional fields specification](https://samtools.github.io/hts-specs/SAMtags.pdf). The default tag is `km:Z:`, followed by the _k_-mers separated by commas. This meets the specification for a string tag and lowercase characters are reserved for user-defined tags.

For a detailed explanation of the matching statistics, see the [**Log Output**](log.md) section.

## Usage

Basic usage of the `tag` subcommand is as follows:

```text
merkurio tag [OPTIONS] -i <IN_FILE> <-s <KMER_SEQ>...|-f <KMER_FILE>>
```

For more examples, see the [**Example Usage**](tag-example.md) section.

## Detailed Description of the Available Options

You can display a help message with `-h` or `--help`. 

### Input/output parameters: 

| Short flag | Long flag     | Description                                                                                                                                                                                                                                                |
| ---------- | ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`       | `--in-file`   | `<Path to input file (SAM/BAM)>` Supports `.bam` and `.sam` files, automatically recognizing the file type based on the extension.                                                                                                                         |
| `-s`       | `--kmer-seq`  | `<Query sequences to search for>...` Multiple sequences can be provided as arguments, separated by spaces. Duplicate sequences are ignored.                                                                                                                |
| `-f`       | `--kmer-file` | `<Path to a file containing query sequences>` Can be in FASTA format or plain text, with empty lines and lines preceded by a `#` being ignored.                                                                                                            |
| `-o`       | `--out-file`  | `<Output file path>` If not provided, output is written to stdout (i.e., the terminal). The extension of the output file path determines the file type (SAM/BAM). If none is provided, the input file type will be used.                                   |
| `-l`       | `--out-log`   | Set this flag without any arguments to write matching statistics to stdout, or write to file if a path to the output file is passed as an argument to this option. For an explanation of the matching statistics, see the [section below](extract-log.md). |
| `-j`       | `--json-log`  | Set this flag without any arguments to write matching statistics in JSON format to stdout, or provide a file path to write JSON log to a file. If both `-l` and `-j` are set without arguments, it will return an error.                                   |
| `-S`       | `--suppress-output`    | Set this flag to suppress the output of matching records. Only the matching statistics are printed (either use `-l` or `-j` for plain text or JSON logging, respectively).                                                                                  |
| `-m`       | `--filter-matching`  | Set this flag to only output records that contain at least one of the _k_-mers. If no _k_-mers are found, the record is not output.                                                                                                                                                                                                                                               |

### Search parameters: 

| Short flag | Long flag              | Description                                                                                                                                                                                                                                                 |
| ---------- | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-r`       | `--reverse-complement` | Set this flag to also search for the reverse complements of the nucleotide query sequences (i.e., a sequence is reversed, and C/G and A/T are swapped). Duplicate sequences are not added to the list. Works for IUPAC codes, everything else is unchanged. |
| `-c`       | `--canonical`          | Set this flag to only search for the canonical forms of k-mers (i.e., the lexicographically first between a sequence and its reverse complement).                                                                                                           |
| `-v`       | `--invert-match`       | Set this flag to invert the matching behavior: instead of keeping records that contain matching _k_-mers, keep only records that do not match any of the _k_-mers. This flag cannot be used together with `-m` (--filter-matching).                                                                                                    |
| `-I`       | `--case-insensitive`   | Set this flag to use case-insensitive matching. Always uses the Aho-Corasick algorithm.                                                                                                                                                                     |
| `-L`       | `--lowercase`          | Set this flag to convert all input sequences to lowercase.                                                                                                                                                                                                  |
| `-U`       | `--uppercase`          | Set this flag to convert all input sequences to uppercase.                                                                                                                                                                                                  |

### Special parameters: 

| Short flag | Long flag        | Description                                                                                                                                                                                                                                                                                                                                                                       |
| ---------- | ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-t`       | `--tag`              | `<Tag to use>` The tag must be exactly two characters long. The default is `km`. Consider the [SAM specifications](https://samtools.github.io/hts-specs/SAMtags.pdf) when choosing a tag to avoid conflicts. If a tag is already present, the new values are appended to existing ones.                                                                                           |
| `-p`       | `--threads`          | `<Number of threads>` The number of parallel threads used when processing BAM files. Default is 1. Has no effect on SAM files.                                                                                                                                                                                                                                                    |
| `-a`       | `--aho-corasick` | Set this flag to use the Aho-Corasick algorithm. It is best used when searching for lots of patterns at once or when searching for _k_-mers with more than 64 characters, which are too long for the BNDMq algorithm.                                                                                                                                                |
| `-q`       | `--q-size`       | `<Size of q-grams for BNDMq>` Use the Backwards Nondeterministic DAWG Matching algorithm tuned with _q_-grams ([BNDMq](https://doi.org/10.1137/1.9781611972894.3)) with this value for the size of _q_. The optimal value of _q_ depends on the size of the query and type of text searched and is usually around 3-5. It must not be larger than the length of the pattern. |
