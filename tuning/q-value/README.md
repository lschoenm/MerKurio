# BNDMq q-value tuning

This directory contains a matcher-only experiment for selecting the BNDMq
q-gram size as a function of pattern length `k`. It deliberately excludes
FASTQ parsing, threading, logging, and output so those costs cannot obscure
the matcher parameter.

The runner compiles the production `BNDMq` source files directly. It tests:

- every `k` from 4 through 64;
- `q = 1..min(k, 12)`;
- 256 distinct, independently generated patterns for every `k`;
- 150-base deterministic DNA sequences;
- 512 sequences per pattern, with the pattern inserted into exactly half;
- both first-match and all-match searches.

Each timed run searches the complete 256-pattern bank, making the reported
ns/base value an aggregate over all patterns. Matcher construction, corpus
generation, and correctness checks against a naive matcher are outside the
timed region. There are fifteen measurement rounds. Both q-value order and pattern
order are freshly shuffled in every round, while every q-value still sees the
same pattern order within that round.

## Run

From this directory:

```bash
./run.sh
```

The default 256-pattern, 15-round sweep took about nine minutes on the machine
described in `results/metadata.txt`.

Runner parameters can be overridden, for example:

```bash
./run.sh --patterns-per-k 16 --sequences 128 --runs 3 --target-ms 5
```

Use `./run.sh --help` for all options. The checked-in results should be
generated with the defaults.

## Results

The `results` directory contains:

- `q_sweep.csv`: raw timing samples;
- `q_summary.csv`: median, mean, spread, and relative performance for every
  `(k, q, mode)`;
- `q_best_by_mode.csv`: independently optimal q-values for first-match and
  all-match searches;
- `q_curves.svg`: two panels with twelve curves, one curve per q-value;
- `q_optimal.svg`: optimal q by `k`, with the percentage advantage over the
  second-best q and a shaded practical-tie region below 1%;
- `metadata.txt`: configuration, compiler, platform, and source revision.

These are internal tuning measurements, not portable performance claims.
Before changing `tune_q_value`, smooth noisy one-base transitions into simple
ranges and validate the proposed ranges on another CPU architecture.
