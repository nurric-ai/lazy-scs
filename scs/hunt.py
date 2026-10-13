"""Search for instances exercising the rare branches of the connection procedure."""
import random, sys
from collections import Counter
from scs2 import two_approx

RARE = {"requests", "host_requests", "individual_link", "below_baseline", "self_cycle", "upper_round_trip"}

def gen(rng):
    out = []
    for _ in range(rng.randint(2, 9)):
        small = rng.choice(["a", "ab", "aab", "abc", "abb", "aabb"])
        tail = rng.choice(["", "c", "b", "cc", "bc", "d", "cd"])
        prim = small * rng.randint(1, 6) + tail
        if rng.random() < 0.3:
            prim = rng.choice(["a", "ab", "aab", "abc"])
        L = rng.randint(2, 40)
        off = rng.randrange(len(prim))
        s = (prim * (L // len(prim) + 3))[off: off + L]
        out.append(s)
    return out

N = int(sys.argv[1]); rng = random.Random(int(sys.argv[2]))
hits, examples = Counter(), {}
for i in range(N):
    inst = gen(rng)
    r = two_approx(inst)
    assert len(r.superstring) <= 2 * r.W and all(s in r.superstring for s in inst)
    for k in RARE:
        if r.stats.get(k):
            hits[k] += 1
            if k not in examples or sum(map(len, inst)) < sum(map(len, examples[k])):
                examples[k] = inst
print(N, "instances OK; rare-branch hits:", dict(hits))
for k, v in examples.items():
    print(k, v)
