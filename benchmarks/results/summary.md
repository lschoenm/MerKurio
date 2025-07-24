# Benchmark Results Summary

## Contents

- [100x31mers-results](#fastq-paired-100x31mers-results)
- [1x31mers-results](#fastq-paired-1x31mers-results)
- [100x31mers-results](#fastq-100x31mers-results)
- [1x31mers-results](#fastq-1x31mers-results)

## fastq-paired: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 1.918 | 0.030 | 1.835 | 1.997 | 1.86x |
| cookiecutter | 8.233 | 0.119 | 7.948 | 8.491 | 7.97x |
| merkurio | 1.033 | 0.086 | 0.919 | 1.172 | 1.00x |

## fastq-paired: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 1.856 | 0.027 | 1.805 | 1.938 | 2.71x |
| cookiecutter | 3.847 | 0.074 | 3.659 | 4.078 | 5.62x |
| merkurio | 0.684 | 0.073 | 0.597 | 0.854 | 1.00x |

## fastq: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 2.609 | 0.058 | 2.511 | 2.792 | 5.08x |
| grep | 1.045 | 0.044 | 0.981 | 1.115 | 2.03x |
| cookiecutter | 2.808 | 0.052 | 2.699 | 2.913 | 5.47x |
| seqkit | 5.531 | 0.156 | 5.345 | 5.932 | 10.77x |
| back_to_sequences | 1.460 | 0.061 | 1.355 | 1.634 | 2.84x |
| merkurio | 0.513 | 0.041 | 0.455 | 0.582 | 1.00x |

## fastq: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 0.430 | 0.043 | 0.373 | 0.519 | 1.19x |
| grep | 0.383 | 0.041 | 0.334 | 0.451 | 1.06x |
| cookiecutter | 1.729 | 0.030 | 1.654 | 1.812 | 4.81x |
| seqkit | 0.891 | 0.055 | 0.791 | 1.070 | 2.48x |
| back_to_sequences | 1.325 | 0.028 | 1.276 | 1.422 | 3.68x |
| merkurio | 0.360 | 0.043 | 0.305 | 0.440 | 1.00x |

