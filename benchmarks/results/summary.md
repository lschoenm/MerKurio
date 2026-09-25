# Benchmark Results Summary

## Contents

- [1000000x31mers-results](#fastq-paired-1000000x31mers-results)
- [100x31mers-results](#fastq-paired-100x31mers-results)
- [1x31mers-results](#fastq-paired-1x31mers-results)
- [1000000x31mers-results](#fastq-1000000x31mers-results)
- [100x31mers-results](#fastq-100x31mers-results)
- [1x31mers-results](#fastq-1x31mers-results)

## fastq-paired: 1000000x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 9.801 | 0.183 | 9.341 | 9.964 | 1.66x |
| merkurio | 5.894 | 0.113 | 5.732 | 6.072 | 1.00x |

## fastq-paired: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 1.851 | 0.055 | 1.719 | 1.976 | 1.39x |
| cookiecutter | 5.958 | 0.080 | 5.790 | 6.177 | 4.47x |
| merkurio | 1.332 | 0.057 | 1.224 | 1.413 | 1.00x |

## fastq-paired: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 1.798 | 0.052 | 1.706 | 1.919 | 2.21x |
| cookiecutter | 3.598 | 0.079 | 3.378 | 3.837 | 4.43x |
| merkurio | 0.812 | 0.059 | 0.736 | 0.947 | 1.00x |

## fastq: 1000000x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| grep | 15.641 | 1.162 | 14.258 | 17.154 | 5.12x |
| cookiecutter | 199.384 | 1.147 | 198.204 | 201.343 | 65.32x |
| back_to_sequences | 4.696 | 0.171 | 4.463 | 5.031 | 1.54x |
| merkurio | 3.052 | 0.096 | 2.960 | 3.234 | 1.00x |

## fastq: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 2.879 | 0.090 | 2.738 | 3.086 | 4.67x |
| grep | 1.200 | 0.022 | 1.158 | 1.251 | 1.94x |
| cookiecutter | 2.579 | 0.054 | 2.461 | 2.720 | 4.18x |
| seqkit | 5.953 | 0.051 | 5.775 | 6.100 | 9.65x |
| back_to_sequences | 1.598 | 0.086 | 1.429 | 1.831 | 2.59x |
| merkurio | 0.617 | 0.034 | 0.567 | 0.678 | 1.00x |

## fastq: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 0.490 | 0.039 | 0.432 | 0.571 | 1.14x |
| grep | 0.480 | 0.041 | 0.388 | 0.548 | 1.11x |
| cookiecutter | 1.844 | 0.056 | 1.739 | 2.068 | 4.28x |
| seqkit | 0.949 | 0.046 | 0.846 | 1.050 | 2.21x |
| back_to_sequences | 1.405 | 0.041 | 1.334 | 1.486 | 3.26x |
| merkurio | 0.431 | 0.055 | 0.335 | 0.516 | 1.00x |

