"""Compressed + lazy: valid superstring, len <= 2W' and W' <= OPT on every instance."""
import random, sys, time
from collections import Counter
from scs2 import reduce_inputs
from baselines import exact
from pipeline_c import two_approx_c_lazy
from stage1 import sources
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 21)
n, worst, rounds, t0 = 0, 0, Counter(), time.time()
for src, inst in sources(rng):
    req = reduce_inputs(inst)
    if not req: continue
    r = two_approx_c_lazy(inst)
    assert all(s in r.superstring for s in inst) and len(r.superstring) <= 2 * r.W
    if len(req) <= 12:
        opt = len(exact(inst)); assert r.W <= opt, (inst, r.W, opt)
        worst = max(worst, len(r.superstring) / opt)
    rounds[r.stats["lazy_rounds"]] += 1; n += 1
print(f"{n} instances OK in {time.time()-t0:.0f}s; worst len/OPT {worst:.3f}; rounds histogram {dict(sorted(rounds.items()))}")
