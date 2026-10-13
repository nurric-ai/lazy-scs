# lazy-scs

Fast, memory-lean **factor-2 approximation for the shortest common superstring** (SCS), with a
**certified lower bound** on the optimum. It implements the first polynomial-time factor-2
algorithm ([openai/math result 128](https://github.com/openai/math)) using lazily enabled periodic
rules and suffix-automaton counts, so memory grows linearly with the input (about 115 bytes per
input character), not quadratically.

```bash
pip install lazy-scs
```

```python
import lazy_scs

lazy_scs.shortest_superstring(["GATTACA", "TACAGG", "AGGCT"])
# 'GATTACAGGCT'

r = lazy_scs.solve(reads)           # list of str or bytes
r.superstring                       # contains every read
r.lower_bound                       # certified: lower_bound <= OPT
r.ratio                             # len / lower_bound; 1.0 = provably optimal
```

The output always satisfies `len(superstring) <= 2 * lower_bound <= 2 * OPT`. On large, realistic
inputs it is usually within a fraction of a percent of the bound; on tiny adversarial inputs it can
be further from optimal (still within the factor 2). `solve(..., verify=True)`
re-justifies every count of the lower bound with an independent checker.

Command line (one string per line, or FASTA):

```bash
lazy-scs reads.txt -o superstring.txt --verify
```

On error-free 150 bp reads at 10x coverage of the *E. coli* K-12 genome (46.4 million input
characters) it runs in about two minutes on a laptop and returns a superstring within 0.021% of the
certified lower bound.

Paper: *Lazy Periodic Rules: A Scalable Implementation of the Factor-2 Approximation for Shortest
Common Superstring*, Krushna Dere and Svyatoblood, independent researchers, October 2026 —
[paper.pdf](https://github.com/nurric-ai/lazy-scs/blob/main/paper/paper.pdf). Source, the Lean 4
formalisation of the lower-bound argument, and benchmarks: <https://github.com/nurric-ai/lazy-scs>.

MIT licence.
