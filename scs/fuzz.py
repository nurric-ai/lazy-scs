"""Fuzz the 2-approximation against the exact optimum.

Checks on every instance: valid superstring, W <= OPT, len <= 2W, and all the
paper's internal invariants (scs2.CHECK). Also tallies which branches of the
connection procedure ran, so we know the hard cases were exercised.

    python fuzz.py [instances] [seed]
"""
import random
import sys
import time
from collections import Counter

from scs2 import two_approx
from baselines import exact


def gen(rng):
    kind = rng.choice(["random", "periodic", "periodic", "mixed", "runs", "nested", "nested"])
    if kind == "nested":
        # periodic texts whose primitive word itself contains a smaller-period run,
        # e.g. (ab)^3 c repeated: forces long agreements between different periods
        out = []
        for _ in range(rng.randint(2, 7)):
            small = rng.choice(["a", "ab", "ba", "aab", "abb"])
            prim = small * rng.randint(1, 4) + rng.choice(["", "b", "c", "ab", "bc"])
            if rng.random() < 0.5:
                prim = small * rng.randint(2, 5)
            L = rng.randint(2, 20)
            off = rng.randrange(len(prim))
            out.append((prim * (L // len(prim) + 3))[off: off + L])
        return out
    n = rng.randint(1, 8)
    sigma = "ab" if rng.random() < 0.6 else "abc"
    out = []
    for _ in range(n):
        L = rng.randint(1, 9)
        if kind == "random":
            s = "".join(rng.choice(sigma) for _ in range(L))
        elif kind == "periodic" or (kind == "mixed" and rng.random() < 0.5):
            P = "".join(rng.choice(sigma) for _ in range(rng.randint(1, 3)))
            off = rng.randrange(len(P))
            s = (P * 12)[off: off + L]
            if rng.random() < 0.3:                      # periodic with a defect
                i = rng.randrange(L)
                s = s[:i] + rng.choice(sigma) + s[i + 1:]
        elif kind == "runs":
            s = "".join(rng.choice(sigma) * rng.randint(1, 4) for _ in range(rng.randint(1, 3)))
        else:
            s = "".join(rng.choice(sigma) for _ in range(L))
        out.append(s)
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    branches, worst, t0 = Counter(), (0, None), time.time()
    for it in range(N):
        inst = gen(rng)
        try:
            r = two_approx(inst)
        except Exception as e:
            print("FAIL on", inst, "->", repr(e))
            raise
        opt = len(exact(inst))
        assert r.W <= opt, ("W > OPT", inst, r.W, opt)
        assert len(r.superstring) <= 2 * r.W
        assert all(s in r.superstring for s in inst)
        for k, v in r.stats.items():
            if k not in ("groups", "layers", "added_cost") and v:
                branches[k] += 1
        ratio = len(r.superstring) / opt if opt else 1
        if ratio > worst[0]:
            worst = (ratio, inst)
    print(f"{N} instances OK in {time.time() - t0:.1f}s; worst len/OPT = {worst[0]:.3f} on {worst[1]}")
    print("instances exercising each branch:", dict(branches))


if __name__ == "__main__":
    main()
