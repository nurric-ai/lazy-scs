
### random, 4 letters  (300 instances, ratio = length / OPT)

| Algorithm | Mean | Worst | Optimal |
|---|---:|---:|---:|
| greedy | 1.005 | 1.074 | 86.7% |
| 2apx | 1.252 | 1.688 | 2.3% |
| 2apx+trim | 1.103 | 1.500 | 34.3% |
| hybrid | 1.004 | 1.074 | 90.0% |

### random, 2 letters  (300 instances, ratio = length / OPT)

| Algorithm | Mean | Worst | Optimal |
|---|---:|---:|---:|
| greedy | 1.012 | 1.136 | 81.3% |
| 2apx | 1.265 | 1.867 | 10.0% |
| 2apx+trim | 1.082 | 1.462 | 44.7% |
| hybrid | 1.007 | 1.133 | 87.7% |

### genome reads  (300 instances, ratio = length / OPT)

| Algorithm | Mean | Worst | Optimal |
|---|---:|---:|---:|
| greedy | 1.000 | 1.026 | 99.7% |
| 2apx | 1.262 | 1.643 | 3.7% |
| 2apx+trim | 1.054 | 1.552 | 73.3% |
| hybrid | 1.000 | 1.026 | 99.7% |

### periodic  (300 instances, ratio = length / OPT)

| Algorithm | Mean | Worst | Optimal |
|---|---:|---:|---:|
| greedy | 1.001 | 1.091 | 97.3% |
| 2apx | 1.441 | 1.939 | 5.0% |
| 2apx+trim | 1.106 | 1.692 | 43.7% |
| hybrid | 1.001 | 1.091 | 97.7% |

### Greedy trap {c(ab)^k, (ba)^k, (ab)^k c}  (lengths)

| k | OPT | greedy | 2apx | 2apx+trim |
|---:|---:|---:|---:|---:|
| 2 | 8 | 10 | 10 | 10 |
| 4 | 12 | 18 | 14 | 14 |
| 8 | 20 | 34 | 22 | 22 |
| 16 | 36 | 66 | 38 | 38 |
| 32 | 68 | 130 | 70 | 70 |

### genome reads, n=300, 40bp reads, 1000bp genome  (5 instances, ratio = length / W, W <= OPT)

| Algorithm | Mean | Worst |
|---|---:|---:|
| greedy | 1.000 | 1.001 |
| 2apx | 1.036 | 1.040 |
| 2apx+trim | 1.000 | 1.001 |
| hybrid | 1.000 | 1.001 |

median time: greedy 9108 ms, 2apx 1543 ms

### random 4 letters, n=200  (5 instances, ratio = length / W, W <= OPT)

| Algorithm | Mean | Worst |
|---|---:|---:|
| greedy | 1.000 | 1.002 |
| 2apx | 1.041 | 1.058 |
| 2apx+trim | 1.026 | 1.046 |
| hybrid | 1.000 | 1.000 |

median time: greedy 3810 ms, 2apx 631 ms
