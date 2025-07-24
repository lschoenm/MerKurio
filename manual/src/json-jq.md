# Useful `jq` Commands

[`jq` is a](https://jqlang.org/) lightweight and flexible command-line JSON processing language. 

A work-in-progress list of `jq` commands to analyse the JSON logs of MerKurio on the command line. Replace `log.json` with the path to your JSON file. 

## Quality metrics

For paired-end reads:

```bash
jq -r '. | "Number of reads with matches: \(.matching_statistics.number_of_distinct_records_with_a_hit)
Number of extracted reads: \(.paired_end_reads_statistics.number_of_extracted_records)
Total number of reads: \(.matching_statistics.number_of_records_searched)
Match rate: \((.matching_statistics.number_of_distinct_records_with_a_hit/.matching_statistics.number_of_records_searched)*100*1000|round/1000)/100% of reads contained target k-mers
Extraction rate: \((.paired_end_reads_statistics.number_of_extracted_records/.matching_statistics.number_of_records_searched)*100*1000|round/1000)/100% of reads were extracted"' log.json
```

## Summary reports

For single files:

```bash
jq -r '
"--- MerKurio Extraction Summary ---",
"Run: " + .meta_information.timestamp,
"Algorithm: " + .meta_information."search-algorithm",
"",
"SEARCH RESULTS:",
"• Patterns searched: " + (.matching_statistics.number_of_patterns_searched | tostring),
"• Patterns found: " + (.matching_statistics.number_of_patterns_found | tostring) + " (" + ((.matching_statistics.number_of_patterns_found/.matching_statistics.number_of_patterns_searched)*100*100 | round/100 | tostring) + "%)",
"• Total matches: " + (.matching_statistics.number_of_matches | tostring),
"• Reads with hits: " + (.matching_statistics.number_of_distinct_records_with_a_hit | tostring) + "/" + (.matching_statistics.number_of_records_searched | tostring),
""
' log.json
```

For paired-end reads:

```bash
jq -r '
"--- MerKurio Extraction Summary ---",
"Run: " + .meta_information.timestamp,
"Algorithm: " + .meta_information."search-algorithm",
"",
"SEARCH RESULTS:",
"• Patterns searched: " + (.matching_statistics.number_of_patterns_searched | tostring),
"• Patterns found: " + (.matching_statistics.number_of_patterns_found | tostring) + " (" + ((.matching_statistics.number_of_patterns_found/.matching_statistics.number_of_patterns_searched)*100*100 | round/100 | tostring) + "%)",
"• Total matches: " + (.matching_statistics.number_of_matches | tostring),
"• Reads with hits: " + (.matching_statistics.number_of_distinct_records_with_a_hit | tostring) + "/" + (.matching_statistics.number_of_records_searched | tostring),
"",
"PAIRED-END DETAILS:",
"• R1 hits: " + (.paired_end_reads_statistics.number_of_hits_in_file_1 | tostring) + " in " + (.paired_end_reads_statistics.number_of_distinct_records_with_a_hit_in_file_1 | tostring) + " reads",
"• R2 hits: " + (.paired_end_reads_statistics.number_of_hits_in_file_2 | tostring) + " in " + (.paired_end_reads_statistics.number_of_distinct_records_with_a_hit_in_file_2 | tostring) + " reads",
"• Extracted pairs: " + (.paired_end_reads_statistics.number_of_extracted_records | tostring),
""
' log.json
```

## Data export

_K_-mers with at least one match to a FASTA file, with the match count in the sequence header:

```bash
jq -r '
  .pattern_hit_counts
  | to_entries
  | map(select(.value > 0))
  | to_entries
  | map(">kmer"+((.key + 1)|tostring)+"|count="+(.value.value|tostring)+"\n"+.value.key)
  | .[]
' log.json > matched_kmers.fasta
```

Pattern hit counts to TSV:

```bash
jq -r '.pattern_hit_counts | to_entries | ["Pattern", "Count"], (.[] | [.key, .value]) | @tsv' log.json > kmer_counts.tsv
```

Pattern positions to TSV:

```bash
jq -r '.matching_records | ["Pattern", "Position", "File", "ReadID"], (.[] | [.pattern, .position, .file, .record_id]) | @tsv' log.json > kmer_positions.tsv
```

## Multiple hits of the same patterns

Producing compact output (empty if nothing is found):

```bash
jq -r '
.matching_records
| group_by(.pattern)
| .[]
| . as $pattern_group
| ($pattern_group | group_by(.record_id) | map(select(length > 1))) as $multi_hits
| if ($multi_hits | length) > 0 then
    "   Pattern: " + $pattern_group[0].pattern,
    "   Reads with multiple hits:",
    ($multi_hits[] | "   • " + .[0].record_id + " (" + (length | tostring) + " hits at positions: " + ([.[].position | tonumber] | sort | map(tostring) | join(", ")) + ")")
  else empty end
' log.json
```

Display information line-wise:

```bash
jq -r '
.matching_records
| group_by(.pattern)
| .[]
| . as $pattern_group
| ($pattern_group | group_by(.record_id) | map(select(length > 1))) as $multi_hits
| if ($multi_hits | length) > 0 then
    "",
    "=" * 60,
    "PATTERN: " + $pattern_group[0].pattern,
    "=" * 60,
    ($multi_hits[] |
      "Read: " + .[0].record_id + " (File: " + .[0].file + ")",
      "Hits: " + (length | tostring),
      (. | sort_by(.position | tonumber) | .[] | "  Position " + .position),
      ""
    )
  else empty end
' log.json
```

Less detailed output:

```bash
jq -r '
.matching_records
| group_by(.pattern)
| map({
    pattern: .[0].pattern,
    multi_hit_reads: (. | group_by(.record_id) | map(select(length > 1)) | length),
    total_multi_hits: (. | group_by(.record_id) | map(select(length > 1)) | map(length) | add // 0)
  })
| map(select(.multi_hit_reads > 0))
| if length > 0 then
    .[] | .pattern[0:30] + "... : " + (.multi_hit_reads | tostring) + " read(s) with " + (.total_multi_hits | tostring) + " total multi-hits"
  else
    "No reads found with multiple hits of the same pattern"
  end
' log.json
```
