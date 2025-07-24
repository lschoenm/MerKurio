# JSON Log

The log in [JSON](https://www.json.org/) format is intended to be machine-readable. It contains the same information as the plain text log but in a structured format. The JSON log contains five main sections: `matching_records`, `summary_statistics`, `meta_information`, `paired_end_reads_statistics`, and `pattern_hit_counts`.

The `matching_records` array contains a list of all matches. Each match is stored as an object with the file name, record ID, query sequence, and zero-based position of the match in that record.

The `summary_statistics` object contains the total number of records searched, the total number of characters of sequences searched, the total number of hits, the number of records with at least one hit, and the numbers of searched/matching patterns.

The `meta_information` object contains the command passed to execute MerKurio as an array, the program's name (MerKurio) and version, the timestamp when the log was generated, and the SAM tag in case of the `tag` subcommand. It also stores the names of input files in an object and information about the search mode (inverted matching extracts only non-matching records, case-insensitive search, used algorithm). 

The `paired_end_reads_statistics` object contains information about the number of hits in each file, the number of records with at least one hit for each file, and the total number of extracted records. In paired-end read mode, a match in one read of a pair will extract both of them. The total number of extracted records can thus be higher than the number of distinct records with a hit. The `searching_paired_end_reads` boolean indicates whether paired-end reads were used (i. e., a second read file was passed via the `-2` flag).

The `pattern_hit_counts` field is a dictionary, with an entry for each pattern searched and the number of times it was found.

```json
{
  "matching_records": [
    {
      "file": "reads_1.fastq",
      "pattern": "ATCG",
      "position": "0",
      "record_id": "record 1/1"
    },
    {
      "file": "reads_2.fastq",
      "pattern": "CGAT",
      "position": "147",
      "record_id": "record 1/2"
    },
    {
      "file": "reads_1.fastq",
      "pattern": "ATCG",
      "position": "13",
      "record_id": "record 1/1"
    }
  ],
  "meta_information": {
    "case_insensitive": false,
    "command_line": [
      "merkurio",
      "extract",
      "-i",
      "reads_1.fastq",
      "-2",
      "reads_2.fastq",
      "-s",
      "ATCG",
      "-r",
      "-j",
      "log.json"
    ],
    "input_files": {
      "kmer_file": null,
      "record_file_1": "reads_1.fastq",
      "record_file_2": "reads_2.fastq"
    },
    "inverted_matching": false,
    "program": "merkurio",
    "search_algorithm": "BNDMq",
    "subcommand": "extract",
    "timestamp": "2024-09-23T00:13:17+02:00[Europe/Vienna]",
    "version": "1.0.0"
  },
  "paired_end_reads_statistics": {
    "number_of_distinct_records_with_a_hit_in_file_1": 2,
    "number_of_distinct_records_with_a_hit_in_file_2": 1,
    "number_of_extracted_records": 4,
    "number_of_hits_in_file_1": 2,
    "number_of_hits_in_file_2": 1,
    "searching_paired_end_reads": true
  },
  "pattern_hit_counts": {
    "ATCG": 2,
    "CGAT": 1
  },
    "summary_statistics": {
    "number_of_characters_searched": 150000,
    "number_of_distinct_records_with_a_hit": 3,
    "number_of_matches": 3,
    "number_of_patterns_found": 2,
    "number_of_patterns_searched": 2,
    "number_of_records_searched": 10000
  }
}
```

This JSON output can easily be read by other programming languages [like Python](./json-python.md), or processed with other command-line tools [such as `jq`](./json-jq.md). For tutorials on both, see the next sections. 