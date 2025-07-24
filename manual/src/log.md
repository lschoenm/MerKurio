# Explanation of the Log Matching Statistics

Recording detailed match statistics is available for both the `extract` and `tag` subcommands. The logs of the two subcommands are mostly identical, with the exception of the `paired_end_reads_statistics` section, which is only present in the `extract` subcommand log after processing paired-end reads, and an additional field `tag` in the [JSON](https://www.json.org/) log of the `tag` subcommand. The examples on this page are based on the `extract` subcommand.

There are two types of log output format available: plain text and JSON. An example of each is provided below. Except for minor differences in formatting, the information contained in both log formats is the same. Plain text and JSON logs are written to files, provided with the `-l` and `-j` arguments, respectively. If only the flag are present, the log will be written to the terminal (stdout). 

When processing paired-end reads, additional statistics are collected. If a single file is processed, the number of extracted records equals the number of distinct records with a hit. In paired-end reads, however, a match in one read of a pair will extract both of them.

## Plain Text Log

The log in plain text format is meant to be human-readable. Lines starting with a `#` contain meta-information and matching statistics. The log also contains a tab-delimited table with an entry for each pattern found in a record, listing the file name, record ID, query sequence and zero-based position of the match. These lines are not prefixed with a `#` to facilitate inspection and further analysis, e.g. using `grep -v '^#'` to extract a tab-delimited table.

After the table, patterns and their number of occurences are listed. This is followed by summary statistics. 

```text
#MerKurio extract log
#2024-09-23T00:13:17+02:00[Europe/Vienna]
#Running MerKurio version 1.0.0
#Command line: merkurio extract -i reads_1.fastq -2 reads_2.fastq -s ATCG -l -r
#Searching for 2 patterns
#
#File           Record       Pattern  Position (zero-based)
reads_1.fastq   record 1/1   ATCG     0
reads_2.fastq   record 1/2   CGAT     147
reads_1.fastq   record 57/1  ATCG     13
#
#Number of patterns found: 2/2 (100.00 %)
#Pattern     Count
#ATCG        2
#CGAT        1
#
#Total number of records searched: 10000
#Total number of characters searched: 150000
#Total number of hits: 3
#Number of distinct records with a hit: 3
```

When processing paired-end reads, the log will contain an additional block of statistics:

```text
#
#Total number of hits in file 1: 2
#Total number of hits in file 2: 1
#Number of distinct records with a hit in file 1: 2
#Number of distinct records with a hit in file 2: 1
#Total number of extracted records: 4
```
