"""Benchmark: paper's 2-approximation vs greedy (and a hybrid of both).

Algorithms
  greedy       classic max-overlap merge (no proven factor below 3.5)
  2apx         openai/math #128 exactly as in the paper (Euler tour output)
  2apx+trim    same, then drop redundant leading/trailing letters
  hybrid       shorter of greedy and 2apx+trim  -> keeps the factor-2 guarantee

    python bench.py
"""
import random
import statistics as st
import time

from scs2 import two_approx
from baselines import exact, greedy, trim


def fam_random(rng, n, L, sigma):
    return ["".join(rng.choice(sigma) for _ in range(rng.randint(max(1, L // 2), L))) for _ in range(n)]


def fam_reads(rng, n, read_len, genome_len):
    g = "".join(rng.choice("ACGT") for _ in range(genome_len))
    out = []
    for _ in range(n):
        i = rng.randrange(genome_len - read_len + 1)
        out.append(g[i:i + read_len])
    return out


def fam_periodic(rng, n, L):
    out = []
    for _ in range(n):
        P = "".join(rng.choice("ab") for _ in range(rng.randint(1, 4)))
        off = rng.randrange(len(P))
        out.append((P * (L // len(P) + 2))[off: off + rng.randint(2, L)])
    return out


def run_algos(inst):
    res, times = {}, {}
    t = time.perf_counter(); g = greedy(inst); times["greedy"] = time.perf_counter() - t
    t = time.perf_counter(); r = two_approx(inst); times["2apx"] = time.perf_counter() - t
    tr = trim(r.superstring, inst)
    res["greedy"], res["2apx"], res["2apx+trim"] = len(g), len(r.superstring), len(tr)
    res["hybrid"] = min(len(g), len(tr))
    for s in inst:
        assert s in g and s in r.superstring and s in tr
    assert len(r.superstring) <= 2 * r.W
    return res, times, r.W


ALGOS = ["greedy", "2apx", "2apx+trim", "hybrid"]


def table_small(name, make, count, rng):
    ratios = {a: [] for a in ALGOS}
    for _ in range(count):
        inst = make(rng)
        res, _, W = run_algos(inst)
        opt = len(exact(inst))
        for a in ALGOS:
            ratios[a].append(res[a] / opt)
    print(f"\n### {name}  ({count} instances, ratio = length / OPT)\n")
    print("| Algorithm | Mean | Worst | Optimal |")
    print("|---|---:|---:|---:|")
    for a in ALGOS:
        v = ratios[a]
        print(f"| {a} | {st.mean(v):.3f} | {max(v):.3f} | {100 * sum(x <= 1 + 1e-9 for x in v) / len(v):.1f}% |")


def table_large(name, make, count, rng):
    ratios, times = {a: [] for a in ALGOS}, {"greedy": [], "2apx": []}
    for _ in range(count):
        inst = make(rng)
        res, tm, W = run_algos(inst)
        for a in ALGOS:
            ratios[a].append(res[a] / W)
        for k in times:
            times[k].append(tm[k])
    print(f"\n### {name}  ({count} instances, ratio = length / W, W <= OPT)\n")
    print("| Algorithm | Mean | Worst |")
    print("|---|---:|---:|")
    for a in ALGOS:
        print(f"| {a} | {st.mean(ratios[a]):.3f} | {max(ratios[a]):.3f} |")
    print(f"\nmedian time: greedy {1e3 * st.median(times['greedy']):.0f} ms, "
          f"2apx {1e3 * st.median(times['2apx']):.0f} ms")


def greedy_trap():
    print("\n### Greedy trap {c(ab)^k, (ba)^k, (ab)^k c}  (lengths)\n")
    print("| k | OPT | greedy | 2apx | 2apx+trim |")
    print("|---:|---:|---:|---:|---:|")
    for k in [2, 4, 8, 16, 32]:
        inst = ["c" + "ab" * k, "ba" * k, "ab" * k + "c"]
        res, _, _ = run_algos(inst)
        print(f"| {k} | {len(exact(inst))} | {res['greedy']} | {res['2apx']} | {res['2apx+trim']} |")


if __name__ == "__main__":
    rng = random.Random(2026)
    table_small("random, 4 letters", lambda r: fam_random(r, r.randint(3, 9), 8, "abcd"), 300, rng)
    table_small("random, 2 letters", lambda r: fam_random(r, r.randint(3, 9), 8, "ab"), 300, rng)
    table_small("genome reads", lambda r: fam_reads(r, r.randint(3, 9), 12, 40), 300, rng)
    table_small("periodic", lambda r: fam_periodic(r, r.randint(3, 8), 12), 300, rng)
    greedy_trap()
    table_large("genome reads, n=300, 40bp reads, 1000bp genome", lambda r: fam_reads(r, 300, 40, 1000), 5, rng)
    table_large("random 4 letters, n=200", lambda r: fam_random(r, 200, 20, "ACGT"), 5, rng)
