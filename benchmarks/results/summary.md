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
| fetch_reads | 2.664 | 1.466 | 1.415 | 4.867 | 1.69x |
| cookiecutter | 50.325 | 0.849 | 48.756 | 51.577 | 31.83x |
| merkurio | 1.581 | 0.105 | 1.488 | 1.862 | 1.00x |

## fastq-paired: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 0.510 | 0.007 | 0.498 | 0.526 | 1.00x |
| cookiecutter | 2.171 | 0.021 | 2.132 | 2.226 | 4.25x |
| merkurio | 0.541 | 0.007 | 0.532 | 0.563 | 1.06x |

## fastq-paired: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| fetch_reads | 0.501 | 0.009 | 0.492 | 0.569 | 5.15x |
| cookiecutter | 1.071 | 0.016 | 1.042 | 1.140 | 11.02x |
| merkurio | 0.097 | 0.002 | 0.095 | 0.107 | 1.00x |

## fastq: 1000000x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| grep | 10.356 | 0.579 | 9.862 | 11.676 | 4.14x |
| cookiecutter | 18.203 | 0.318 | 17.712 | 18.507 | 7.28x |
| back_to_sequences | 3.066 | 0.034 | 3.026 | 3.130 | 1.23x |
| merkurio | 2.499 | 0.860 | 0.955 | 3.177 | 1.00x |

## fastq: 100x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 2.341 | 0.070 | 2.222 | 2.455 | 12.35x |
| grep | 0.745 | 0.029 | 0.707 | 0.800 | 3.93x |
| cookiecutter | 0.983 | 0.039 | 0.928 | 1.046 | 5.19x |
| seqkit | 4.647 | 0.075 | 4.535 | 4.919 | 24.53x |
| back_to_sequences | 1.243 | 0.081 | 1.102 | 1.512 | 6.56x |
| merkurio | 0.189 | 0.004 | 0.183 | 0.201 | 1.00x |

## fastq: 1x31mers-results

| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|---:|
| seqtool | 0.065 | 0.003 | 0.061 | 0.074 | 1.46x |
| grep | 0.084 | 0.006 | 0.077 | 0.105 | 1.87x |
| cookiecutter | 0.518 | 0.007 | 0.509 | 0.536 | 11.57x |
| seqkit | 0.138 | 0.005 | 0.132 | 0.150 | 3.08x |
| back_to_sequences | 1.072 | 0.027 | 1.043 | 1.156 | 23.95x |
| merkurio | 0.045 | 0.004 | 0.040 | 0.059 | 1.00x |

