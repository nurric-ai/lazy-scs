import random, time, tracemalloc, sys
from scs2 import Counts, reduce_inputs
from compressed import CompressedCounts
from bench import fam_reads
rng = random.Random(7)
print("| reads x len | genome | L | ref time | ref mem | compressed time | compressed mem | states | breakpoints | W |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for n, rl, G, do_ref in [(300, 60, 3000, 1), (1000, 100, 10000, 1), (3000, 150, 30000, 0), (10000, 150, 100000, 0), (30000, 150, 300000, 0)]:
    inst = fam_reads(rng, n, rl, G)
    req = reduce_inputs(inst)
    L = sum(map(len, req))
    row = []
    if do_ref:
        tracemalloc.start(); t = time.time()
        ref = Counts(req, rule_filter=lambda *a: False)
        row += [f"{time.time()-t:.1f}s", f"{tracemalloc.get_traced_memory()[1]/2**20:.0f} MB"]
        tracemalloc.stop(); del ref
    else:
        row += ["-", "-"]
    tracemalloc.start(); t = time.time()
    cc = CompressedCounts(req)
    row += [f"{time.time()-t:.1f}s", f"{tracemalloc.get_traced_memory()[1]/2**20:.0f} MB"]
    tracemalloc.stop()
    print(f"| {len(req)} x {rl} | {G} | {L} | " + " | ".join(row) + f" | {cc.states} | {cc.breakpoints} | {cc.W} |", flush=True)
