# Introduction

[MerKurio](https://github.com/lschoenm/MerKurio) was developed as command line tool and provides two main functionalities through dedicated subcommands.

The `extract` subcommand identifies and extracts sequence records from fasta or fastq files based on query _k_-mers, with full support for paired-end reads where matching in one read extracts both mates.

The `tag` subcommand annotates sequences in SAM or BAM alignment formats according to the _k_-mers that they contain following [the SAM Optional Fields Specification](https://samtools.github.io/hts-specs/SAMtags.pdf) with user-defined labels (default "km").

It was developed to sinmplify downstream analysis of selected _k_-mers by tracing them back to their sequences in the original data.

For installation instructions, go to the [installation](./installation.md) section.

To test MerKurio on a minimal working example, visit the [minimal example](./minimal.md) section.

For an example workflow, visit the [practical tutorial](./practical-tutorial.md) section.

You can use the **search bar** on top to find specific information, or use the navigation on the left to browse the manual.

If you have found a bug or have a feature request, please open an issue on the [GitHub repository](https://github.com/lschoenm/MerKurio/issues) or check the [discussions on GitHub](https://github.com/lschoenm/MerKurio/discussions).

## License

MerKurio is licensed under the MIT license. See the [LICENSE](https://github.com/lschoenm/MerKurio/blob/master/LICENSE).
