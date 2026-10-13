"""Stage 1b: does the lazy-rule variant hold up, and how many rules does it need?

For every instance: lazy output must be a superstring with len <= 2W' (hence
<= 2 OPT, since W' <= OPT). We record how many periodic rules were added, how many
rounds it took, and compare W' and output length with the full rule set.

    python stage1_lazy.py [seed]
"""
import random
import statistics as st
import sys
import time
from collections import Counter, defaultdict

from scs2 import Counts, reduce_inputs, two_approx, two_approx_lazy
from baselines import exact
from stage1 import sources


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
    by_src = defaultdict(lambda: defaultdict(list))
    fails = []
    t0 = time.time()
    for src, inst in sources(rng):
        req = reduce_inputs(inst)
        if not req:
            continue
        try:
            lz = two_approx_lazy(inst)
        except Exception as e:
            fails.append((src, inst, repr(e)))
            continue
        assert all(s in lz.superstring for s in inst) and len(lz.superstring) <= 2 * lz.W
        full = two_approx(inst)
        C = Counts(req)
        d = by_src[src]
        d["n_rules_full"].append(C.n_rules)
        d["lazy_rules"].append(lz.stats["lazy_rules"])
        d["lazy_rounds"].append(lz.stats["lazy_rounds"])
        d["W_ratio"].append(lz.W / full.W)
        d["len_ratio"].append(len(lz.superstring) / len(full.superstring))
        if len(req) <= 12:
            opt = len(exact(inst))
            d["lazy_len_opt"].append(len(lz.superstring) / opt)
            d["full_len_opt"].append(len(full.superstring) / opt)
    print(f"done in {time.time() - t0:.0f}s; failures: {len(fails)}")
    for f in fails[:5]:
        print("FAIL", f)
    print("\n| source | n | rules (full set) | rules added lazily: mean / max | instances needing 0 | rounds max | W'/W | lazy len / full len | lazy len/OPT (mean, worst) | full len/OPT (mean, worst) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for src, d in by_src.items():
        n = len(d["lazy_rules"])
        z = sum(1 for x in d["lazy_rules"] if x == 0)
        lo = d.get("lazy_len_opt") or [0]
        fo = d.get("full_len_opt") or [0]
        print(f"| {src} | {n} | {st.mean(d['n_rules_full']):.0f} | {st.mean(d['lazy_rules']):.2f} / {max(d['lazy_rules'])} | "
              f"{100 * z / n:.0f}% | {max(d['lazy_rounds'])} | {st.mean(d['W_ratio']):.3f} | {st.mean(d['len_ratio']):.3f} | "
              f"{st.mean(lo):.3f}, {max(lo):.3f} | {st.mean(fo):.3f}, {max(fo):.3f} |")


if __name__ == "__main__":
    main()
