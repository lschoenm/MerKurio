# Pattern-matching algorithm selection tuning

This directory contains a matcher-only experiment for selecting between
BNDMq, rolling-hash matching, and Aho–Corasick as a function of pattern length
`k` and the number of patterns. It compiles the production matcher source
directly and excludes parsing, threading, logging, and output.

The primary sweep uses deterministic 150-base DNA sequences. Half of the
sequences receive one inserted pattern; naturally occurring additional
matches are retained. Every algorithm sees the same pattern banks and
sequences. The two measured operations are:

- `first`: stop after finding any pattern;
- `all`: enumerate every overlapping occurrence.

Pattern counts are logarithmically spaced up to one million. Counts that
exceed the number of distinct DNA patterns possible for a given `k` are
recorded as invalid rather than run.

## Timeouts and pruning

Every `(k, pattern count, algorithm, mode)` cell runs in a separate child
process. This is necessary because an in-process timer cannot interrupt an
expensive matcher construction or a single search call safely.

Matcher construction is outside the search-throughput timing, but its duration
is recorded separately as `build_ns`. The hard cell timeout covers both
construction and searching.

Two limits apply:

- `--max-sample-ms` rejects a cell when its initial single-iteration search is
  already too slow;
- `--cell-timeout-seconds` lets the parent terminate a stuck or excessively
  expensive worker.

After an algorithm times out at a pattern count, larger counts for the same
`(k, mode)` are marked `pruned`. BNDMq cells above the machine word size are
marked invalid. Timeouts, errors, invalid cells, and pruning decisions are
stored in `results/cell_status.csv`.

## Run

From this directory:

```bash
./run.sh
```

The defaults test:

- `k = 4,6,8,10,12,16,20,24,31,40,48,64,80,96,128`;
- pattern counts
  `1,2,4,8,16,32,64,128,256,1000,10000,100000,1000000`;
- 256 sequences per independent pattern bank;
- up to 64 banks, chosen adaptively to test more pattern compositions at
  small pattern counts;
- 12 measurement rounds, targeting at least 20 ms per timing;
- a 1-second soft sample limit and 30-second hard cell limit.

A small validation run can be launched with:

```bash
./run.sh \
  --k-values 8,31,80 \
  --pattern-counts 1,16,256 \
  --algorithms bndmq,hash,aho_corasick \
  --modes first,all \
  --sequences 32 \
  --runs 2 \
  --target-ms 2 \
  --max-sample-ms 100 \
  --cell-timeout-seconds 5
```

Use `./run.sh --help` for all options.

## Results

- `algorithm_sweep.csv`: raw successful timing samples;
- `cell_status.csv`: completion, timeout, pruning, invalid, and error status;
- `algorithm_summary.csv`: median and spread for each successful algorithm;
- `algorithm_winners.csv`: best and second-best algorithms with their margin;
- `selection_map.svg`: winner heatmaps for first- and all-match search;
- `metadata.txt`: configuration, compiler, platform, and source revision.

The analyzer verifies that successful algorithms produce the same untimed
result checksum for each cell. A mismatch aborts analysis instead of producing
a recommendation.

The currently checked-in result files are the small end-to-end validation
sweep described above, not a production tuning run or recommendation. Replace
them by running the complete default sweep before deriving selection rules.

The winner map is evidence for a later, deliberately simple production
selection rule. It should not be copied directly into a large lookup table.
