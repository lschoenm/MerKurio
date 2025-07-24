pub mod cmd_extract;
pub mod cmd_tag;
pub mod helpers;
pub mod logger;
pub mod pattern_matching;
pub mod pattern_preprocessing;

use anyhow::Result;
use clap::{Parser, Subcommand, crate_authors, crate_version};

use std::env;
use std::str;

#[derive(Parser)]
#[clap(
    author(crate_authors!()),
    version(crate_version!()),
    about("SeqKatcher has two subcommands, 'extract' and 'tag'. The 'extract' subcommand searches for query sequences in FASTA/Q files and extracts records containing the patterns. The 'tag' subcommand filters and tags records in a SAM/BAM file with the presence of query sequences.\n\nThe full documentation can be found here: <https://lschoenm.github.io/MerKurio/>."),
    long_about("SeqKatcher has two subcommands, 'extract' and 'tag'.\nThe 'extract' subcommand searches for query sequences in FASTA/Q files and extracts records containing the patterns. It can also generate a log file with match statistics in plain text or JSON format. Paired-end reads can be processed together.\n\nThe 'tag' subcommand annotates records in a SAM/BAM file with the presence of query sequences by using a SAM optional tag (default 'km'). It can also generate a log with matching statistics.\n\nThe full documentation can be found here: <https://lschoenm.github.io/MerKurio/>.\n\nFor more information, visit the GitHub repository <https://github.com/lschoenm/MerKurio>."),
    arg_required_else_help = true
)]
struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[clap(
        name = "extract",
        about = "Search for query sequences in FASTA/Q files and extract records containing the patterns",
        long_about = "Search for query sequences (k-mers) in FASTA/FASTQ files and extract records containing the patterns. It can also generate a log file with match statistics in plain text or JSON format. Paired-end reads can be processed together.",
        arg_required_else_help = true
    )]
    Extract(cmd_extract::CmdExtract),
    #[clap(
        name = "tag",
        about = "Tag records in a BAM/SAM file with the presence of query sequences",
        long_about = "Tag and filter records in a BAM/SAM file with the presence of query sequences by using a SAM optional tag (default 'km'). Optionally, keep only those records with/without matching k-mers. It can also generate a log file with match statistics in plain text or JSON format.",
        arg_required_else_help = true
    )]
    Tag(cmd_tag::CmdTag),
}

fn main() -> Result<()> {
    // Parse command line arguments
    let args_parsed = Cli::parse();

    // Call subcommand function
    match args_parsed.cmd {
        Commands::Extract(args) => cmd_extract::extract_records(args),
        Commands::Tag(args) => cmd_tag::tag_records(args),
    }
}

//
// ---------------------------------- Tests ----------------------------------
//

#[cfg(test)]
mod tests {
    use super::*;
    use clap::crate_name;

    // Test the Cli parser with common arguments (program not executed)
    #[test]
    fn test_cli_parser_common_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample_1.fasta",
                "-2",
                "tests/data/sample_2.fasta",
                "--kmer-file",
                "tests/data/kmers.txt",
                "-o",
                "tests/data/output.fasta",
                "-l",
                "-r",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    fn test_cli_parser_common_tag() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "tag",
                "-i",
                "tests/data/sample.sam",
                "--kmer-file",
                "tests/data/kmers.txt",
                "-o",
                "tests/data/output.sam",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_group_kmers_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.fasta",
                "--kmer-file",
                "tests/data/kmers.txt",
                "--kmer-seq",
                "ACGT",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_group_kmers_tag() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "tag",
                "-i",
                "tests/data/sample.sam",
                "--kmer-file",
                "tests/data/kmers.txt",
                "--kmer-seq",
                "ACGT",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_group_algorithm_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.fasta",
                "--kmer-file",
                "tests/data/kmers.txt",
                "--aho-corasick",
                "--q-size",
                "5",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    fn test_cli_parser_multiple_seq_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.fasta",
                "--kmer-seq",
                "ACGT",
                "TGCA",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    fn test_cli_parser_multiple_seq_tag() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "tag",
                "-i",
                "tests/data/sample.sam",
                "-o",
                "tests/data/output.sam",
                "--kmer-seq",
                "ACGT",
                "TGCA",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    fn test_cli_parser_suppress_output_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.sam",
                "--kmer-seq",
                "ACGT",
                "-S",
                "-l",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_suppress_output_requires_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.sam",
                "--kmer-seq",
                "ACGT",
                "--suppress-output",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_suppress_output_conflicts_extract() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.sam",
                "-o",
                "tests/data/output.sam",
                "--kmer-seq",
                "ACGT",
                "--suppress-output",
                "-l",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    fn test_cli_parser_suppress_and_json() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.sam",
                "--kmer-seq",
                "ACGT",
                "-S",
                "-j",
                "tests/data/log.json",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }

    #[test]
    #[should_panic]
    fn test_cli_parser_suppress_requires_logging() {
        let args = Cli::try_parse_from(
            vec![
                crate_name!(),
                "extract",
                "-i",
                "tests/data/sample.sam",
                "--kmer-seq",
                "ACGT",
                "-S",
            ]
            .iter(),
        );
        assert!(args.is_ok());
    }
}
