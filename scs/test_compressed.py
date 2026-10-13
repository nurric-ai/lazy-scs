"""Stage 2a check: compressed counts == reference extension-only counts, word by word."""
import random, sys, time
from scs2 import Counts, reduce_inputs
from compressed import CompressedCounts, ev
from stage1 import sources

NORULES = lambda v, p, k: False
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
n = 0
for src, inst in sources(rng):
    req = reduce_inputs(inst)
    if not req:
        continue
    ref = Counts(req, rule_filter=NORULES)
    cc = CompressedCounts(inst)
    assert cc.W == ref.W, (src, inst, cc.W, ref.W)
    up, down = ref.edges()
    S = cc.sam
    for s, ms in ref.m.items():
        if not s:
            continue
        X = S.state_of(s); l = len(s)
        assert cc.m(s) == ms, (src, inst, s, cc.m(s), ms)
        hi = S.len[X]
        u = ms - (ev(cc.f[X], l + 1) if l < hi else cc.children_sum[X])
        d = ms - ev(cc.R[X], l)
        assert (u, d) == (up.get(s, 0), down.get(s, 0)), (inst, s, (u, d), (up.get(s, 0), down.get(s, 0)))
    n += 1
print(n, "instances: counts, W and edge multiplicities identical to the reference")
