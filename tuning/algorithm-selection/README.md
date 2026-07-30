# Pattern-matching algorithm selection tuning

This directory contains an end-to-end matcher experiment for selecting between
BNDMq, rolling-hash matching, and Aho–Corasick as a function of pattern length
`k` and the number of patterns. It compiles the production matcher source
directly. The primary metric includes matcher construction and searching the
configured corpus. Common CLI work—input parsing, threading, logging, and
output—is outside the experiment.

The primary sweep uses deterministic 150-base DNA sequences. Half of the
sequences receive one inserted pattern; naturally occurring additional
matches are retained. Every algorithm sees the same pattern banks and
sequences. The two measured operations are:

- `first`: stop after finding any pattern;
- `all`: enumerate every overlapping occurrence.

Pattern counts are logarithmically spaced up to one million. Counts that
exceed the number of distinct DNA patterns possible for a given `k` are
recorded as invalid rather than run.

Synthetic patterns and sequences are prepared before timing. A timed
end-to-end iteration constructs a fresh matcher for every independent pattern
bank and searches the cell's distributed corpus once. Each cell receives
approximately the same total base budget regardless of bank count. The raw
results retain construction, search, and total durations separately.

Pattern banks are workload replicates, not additional matchers in one
production invocation. The primary reference time therefore combines the
average fixed cost of one matcher with the directly measured search across the
complete cell corpus. All winner and crossover outputs use this reference
time.

## Timeouts and pruning

Every `(k, pattern count, algorithm, mode)` cell runs in a separate child
process. This is necessary because an in-process timer cannot interrupt an
expensive matcher construction or a single search call safely. The hard cell
timeout covers both construction and searching.

Two limits apply:

- `--max-sample-ms` rejects a cell when its initial single-pass matcher build
  and search are already too slow;
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
- approximately 100 million bases total per cell;
- 150 bases per sequence, distributed evenly across the cell's banks;
- up to 64 banks, chosen adaptively to test more pattern compositions at
  small pattern counts;
- five single-pass measurement rounds;
- a 2-second soft sample limit and 60-second hard cell limit.

Thus, the default selection rule is tuned directly for a representative
100-million-base production workload. `--total-bases` can select a different
fixed reference workload.

A small validation run can be launched with:

```bash
./run.sh \
  --k-values 4,8,31,80 \
  --pattern-counts 1,16,256,1000 \
  --algorithms bndmq,hash,aho_corasick \
  --modes first,all \
  --total-bases 153600 \
  --runs 2 \
  --max-sample-ms 100 \
  --cell-timeout-seconds 5
```

Use `./run.sh --help` for all options.

## Refine algorithm transitions

After completing the full sweep, run one adaptive refinement wave on the same
machine:

```bash
./run-refinement.sh
```

The refinement driver reads `algorithm_winners.csv` and finds adjacent
measurements whose winning algorithms differ. It benchmarks all three
algorithms at:

- the geometric midpoint between adjacent pattern counts at the same `k`;
- the arithmetic midpoint between adjacent pattern lengths at the same pattern
  count.

The BNDMq word-size boundary is anchored directly at `k=65` rather than
bisected. Counts or pattern lengths that are already adjacent integers cannot
be refined. The original raw sweep is never modified.

Run the command again to narrow the remaining transition brackets by one more
level, or request multiple levels at once:

```bash
./run-refinement.sh --waves 2
```

To inspect the next wave without running it:

```bash
./run-refinement.sh --dry-run
```

Each cell is stored independently under `results/refinement/parts`, making an
interrupted refinement resumable. Successful and timed-out parts are both
considered complete. The driver rebuilds the aggregate refinement CSVs and
runs the analyzer after every wave.

Refinement timings must be collected on the same machine as the base sweep.
The driver checks the recorded system and refuses to mix machines unless
`--allow-different-system` is explicitly supplied. It also records the base
sweep timestamp and refuses to combine refinement parts with a subsequently
replaced base sweep. Archive or remove `results/refinement` before refining a
new full sweep.

## Results

- `algorithm_sweep.csv`: raw build, search, and total timing samples;
- `cell_status.csv`: completion, timeout, pruning, invalid, and error status;
- `algorithm_summary.csv`: median per-matcher build, search, and reference
  times plus reference timing spread for each successful algorithm;
- `algorithm_winners.csv`: best and second-best algorithms with their margin;
- `selection_map.svg`: categorical winner heatmaps for the original full-sweep
  grid;
- `refined_selection_map.svg`: scatterplots of combined base and refinement
  winners on continuous pattern-length and logarithmic pattern-count axes;
  circles are base measurements and diamonds are refinement measurements.
  Faint background cells use the four original-grid corner winners: unanimous
  cells receive one color, while disagreements split the cell into four
  corner-colored quadrants;
- `crossover_curves.svg`: per-`k` curves showing each algorithm's runtime
  relative to the fastest available algorithm as pattern count increases
  (extreme slowdowns are visibly capped to preserve detail near the crossover);
- `corpus_size_models.csv`: fixed end-to-end matcher cost and search rate used
  to model other corpus sizes;
- `corpus_size_regions.csv`: lower-envelope algorithm regions over searched
  bases per matcher;
- `corpus_size_crossovers.csv`: modeled transitions between optimal algorithms
  with round-resampling bootstrap intervals and transition support;
- `corpus_size_crossovers.svg`: heatmaps of the first corpus-size transition
  for each `k`, pattern count, and search mode;
- `metadata.txt`: configuration, compiler, platform, and source revision.

The refinement-specific raw data are kept in:

- `refinement/manifest.csv`: transition bracket and completion status for each
  selected midpoint;
- `refinement/algorithm_sweep.csv`: aggregated successful refinement timings;
- `refinement/cell_status.csv`: aggregated refinement statuses;
- `refinement/parts/`: resumable per-cell raw results and metadata.

The analyzer verifies that successful algorithms produce the same untimed
result checksum for each cell. A mismatch aborts analysis instead of producing
a recommendation.

The winner tables and refined scatterplot combine the base and refinement
measurements. Runtime crossover curves and all corpus-size crossover outputs
use only the original regular parameter grid. This keeps sparse, targeted
refinement cells from cluttering those plots or changing the corpus-size
models.

The current result files predate the 100-million-base corpus-budget design.
Rerun the complete sweep before deriving selection rules from this version.
Results remain hardware- and workload-dependent and should not be treated as a
portable recommendation.

The winner map is evidence for a later, deliberately simple production
selection rule. It should not be copied directly into a large lookup table.

## Corpus-size crossover model

For each algorithm, the analyzer models the matcher-dependent runtime for `B`
searched bases as:

```text
T(B) = C + R * B
```

`C` is the median fixed cost per matcher, calculated from total time minus
search time. It therefore includes construction, destruction, and fixed timing
overhead. `R` is the median measured search time per base. The reported regions
are the lower envelope of the available algorithm lines.

The crossover table bootstraps complete timing-round observations 1,000 times.
Its interval is conditional on observing the same algorithm transition;
`bootstrap_transition_support_pct` reports how often that transition appeared.
Low support means the identity or existence of the crossover is unstable.

This is an extrapolation from the configured corpus size, match density,
sequence length, and alphabet. Search rate may not remain perfectly constant
across very different corpus sizes or hardware cache regimes. Timed-out and
pruned algorithms cannot appear on the lower envelope, so “throughout” means
throughout the modeled range among successfully measured algorithms. The plot
adds a successful/structurally-eligible count when some eligible algorithms
were unavailable.
