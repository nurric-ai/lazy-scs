# lazy-scs

[![PyPI](https://img.shields.io/pypi/v/lazy-scs)](https://pypi.org/project/lazy-scs/)
[![wheels](https://github.com/nurric-ai/lazy-scs/actions/workflows/wheels.yml/badge.svg)](https://github.com/nurric-ai/lazy-scs/actions/workflows/wheels.yml)
[![Code DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23362570.svg)](https://doi.org/10.5281/zenodo.23362570)
[![Paper DOI](https://img.shields.io/badge/paper-10.5281%2Fzenodo.23348394-blue)](https://doi.org/10.5281/zenodo.23348394)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

A scalable implementation of the polynomial-time **factor-2 approximation for shortest common
superstring** (SCS), [openai/math result 128](https://github.com/openai/math), using **lazy periodic rules**
and **suffix-automaton counts**.

> **Lazy Periodic Rules: A Scalable Implementation of the Factor-2 Approximation for Shortest Common
> Superstring**, Krushna Dere and Svyatoblood, independent researchers, October 2026. [paper/paper.pdf](paper/paper.pdf) ·
> [doi:10.5281/zenodo.23348394](https://doi.org/10.5281/zenodo.23348394)

Given strings `s1 … sn`, SCS asks for a shortest string containing every `si` as a substring. The
factor-2 algorithm builds its lower bound from forced occurrence counts on every distinct substring,
so other implementations use memory that grows quadratically with the input. This one does not:

- **Lazy periodic rules.** Every subset of the paper's periodic rules (2.6) still gives a valid lower
  bound `W' ≤ OPT`, and the construction invokes a periodic rule at exactly one step, which is checked
  at run time. The rules are therefore enabled only for the periodic texts that need them, and the
  output still satisfies `|T| ≤ 2W' ≤ 2·OPT`.
- **Compressed counts.** Without the periodic rules, the counts on each state of a generalised suffix
  automaton form a nonincreasing step function of the word length. Words are (state, length) pairs, so
  substrings are never stored.
- **Machine-checked core.** The lower-bound argument is formalised in Lean 4 with Mathlib
  ([lean/LazySCS.lean](lean/LazySCS.lean)) and depends only on the standard axioms.
- **Independent checker.** `--verify` re-justifies every count from scratch, which is exactly the
  hypothesis of the formalised `certificate` theorem.

## Results

Error-free 150 bp reads at 10× coverage from the *E. coli* K-12 genome (NC_000913.3), 300 s limit
per run, on a laptop:

| genome prefix | input chars | JavaScript impl. | C++ impl. | **lazy-scs** | length / W |
|---|---:|---|---|---|---:|
| 200 bp | 1,950 | 0.12 s | 111 s | <0.01 s | 1.000 |
| 30 kb | 0.3 M | 24 s, 3.4 GB | timeout | 0.52 s, 44 MB | 1.000 |
| 100 kb | 1.0 M | timeout, 10.9 GB | — | 1.9 s, 143 MB | 1.00001 |
| 4.64 Mb (full) | 46.4 M | — | — | 121 s, 6.3 GB | **1.00021** |

`W` is a certified lower bound on the optimum, so the full-genome output is provably within 0.021% of
the shortest possible superstring. Wherever both finish, the output length equals the JavaScript
implementation's exactly. See the paper for correctness tests against exact optima, branch coverage,
and greedy versus the factor-2 algorithm on Shibata's counterexample to the greedy conjecture.

**What this is not.** It is not a genome assembler: SCS merges repeats longer than the read length and
has no model of sequencing errors. On typical inputs, plain greedy produces superstrings that are just
as short. The contribution is scale plus a worst-case guarantee with a per-instance certificate.

## Install

```bash
pip install lazy-scs
```

Prebuilt wheels for Linux, macOS and Windows (Python 3.9+); other platforms build from source and need Rust.

```python
import lazy_scs

lazy_scs.shortest_superstring(["GATTACA", "TACAGG", "AGGCT"])   # 'GATTACAGGCT'

r = lazy_scs.solve(reads, verify=True)   # list of str or bytes
r.superstring, r.lower_bound, r.ratio    # ratio = len / lower_bound; 1.0 = provably optimal
```

The package also installs a `lazy-scs` command (`lazy-scs reads.txt -o superstring.txt --verify`).

## Usage (Rust command line)

```bash
cd scs/scs-rs
cargo build --release
./target/release/scs-rs reads.txt -o superstring.txt          # one string per line, or FASTA
./target/release/scs-rs reads.txt --verify                    # also run the independent checker
./target/release/scs-rs --batch instances.txt --check         # instances separated by '---' lines
```

The program prints JSON with `W` (the lower bound), `length`, `length_over_W`, the number of lazy-rule
rounds, and which branches of the construction ran. On Windows, if the project path is long, build with
`CARGO_TARGET_DIR` set to a short folder.

## Repository layout

| path | contents |
|---|---|
| `pyproject.toml`, `python/` | the `lazy-scs` Python package (PyO3 bindings in `scs/scs-rs/src/python.rs`) and its tests |
| `.github/workflows/wheels.yml` | builds and tests wheels on every push; publishes to PyPI on a `v*` tag |
| `scs/scs-rs/` | Rust implementation: `sam.rs` (automaton, trimming), `counts.rs` (step-function counts, lazy rules), `approx.rs` (layers, connection, Euler tour), `verify.rs` (independent checker) |
| `scs/scs2.py` | Python reference implementation of the paper, one entry per substring, invariants asserted |
| `scs/compressed.py`, `scs/pipeline_c.py` | Python compressed counts and lazy pipeline |
| `scs/baselines.py` | greedy, exact subset DP, trimming |
| `scs/fuzz.py`, `scs/stage1*.py`, `scs/test_*.py` | correctness and differential tests |
| `scs/hunt_rs.py`, `scs/triggers.json`, `scs/trig_rev.json` | rare-branch search and the trigger instances |
| `scs/shibata.py` | Shibata's greedy counterexample, with the greedy run replayed |
| `scs/bench.py`, `scs/h2h.py`, `scs/h2h_results_*.json` | quality benchmark and head-to-head comparison |
| `scs/LAZY_RULES.md` | the lazy-rule argument in prose |
| `lean/` | Lean formalisation, the Kaggle script used to check it, and the compiler output |
| `paper/` | the paper (HTML source and PDF) |
| `misc/` | an unrelated Strassen benchmark from the same session |

## Checking the Lean proof

```bash
lake +leanprover-community/mathlib4:lean-toolchain new lazyscs math
cd lazyscs && lake exe cache get
cp ../lean/LazySCS.lean . && lake env lean LazySCS.lean
```

Expected output: warnings only, and each `#print axioms` lists `propext`, `Classical.choice`,
`Quot.sound`. [lean/lean_check_output.txt](lean/lean_check_output.txt) is the recorded run (Lean
4.35.0-rc4). The facts the paper proves about true occurrence counts enter as hypotheses.

## Limitations

- The claim that the construction invokes a periodic rule only at its Section 5.3 is an informal argument,
  supported by run-time checks of every step.
- Two branches of the construction (layers below a baseline, single-block cycles) never fired in any
  test.
- The Rust implementation and the checker are tested, not verified.

## Disclosure of AI assistance

An AI assistant (referred to here as wigip) wrote the implementation, ran the experiments and drafted
the paper under the authors' direction. The authors are responsible for the content.

## Credits

The algorithm is OpenAI's (openai/math result 128). Shibata's counterexample is arXiv:2609.01365.
Other implementations: [buberlo/veriscs](https://github.com/buberlo/veriscs) (Lean),
[bkarak/shortest-superstring-js](https://github.com/bkarak/shortest-superstring-js) (JavaScript),
[spele1100/scs-2approx](https://github.com/spele1100/scs-2approx) (C++),
[rongou/superstring-rs](https://github.com/rongou/superstring-rs) (Rust, partial).

## Citation

See [CITATION.cff](CITATION.cff). Licence: MIT.
