"""Differential test: Rust binary vs Python compressed counts + exact optimum."""
import random, subprocess, sys, os, time
from scs2 import reduce_inputs
from compressed import CompressedCounts
from baselines import exact
from stage1 import sources

BIN = os.path.join(os.environ["TEMP"], "scsrs", "release", "scs-rs.exe")
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 31)
insts = []
for src, inst in sources(rng):
    # map non-ASCII symbols (Shibata's b_i) to distinct ASCII letters for byte-level Rust
    table = {}
    def asc(s):
        return "".join(c if ord(c) < 128 else table.setdefault(c, "BCDEFGHIJKLMNOPQRSTUVWXYZ"[len(table)]) for c in s)
    inst = [asc(s) for s in inst]
    if reduce_inputs(inst):
        insts.append(inst)
open("batch.txt", "w").write("\n---\n".join("\n".join(i) for i in insts) + "\n")
t = time.time()
out = subprocess.run([BIN, "--batch", "--check", "batch.txt"], capture_output=True, text=True)
rust_secs = time.time() - t
if out.returncode:
    print(out.stderr[-2000:]); sys.exit(1)
lines = out.stdout.splitlines()
bad, worst, rounds = 0, 0, {}
for i, inst in enumerate(insts):
    w_ext, w, ln, rnd = map(int, lines[2 * i].split()[:4])
    T = lines[2 * i + 1]
    assert len(T) == ln <= 2 * w and all(s in T for s in inst)
    py = CompressedCounts(inst).W
    if w_ext != py:
        bad += 1; print("W_ext mismatch", inst, w_ext, py)
    if len(reduce_inputs(inst)) <= 12:
        opt = len(exact(inst)); assert w <= opt, (inst, w, opt); worst = max(worst, ln / opt)
    rounds[rnd] = rounds.get(rnd, 0) + 1
print(f"{len(insts)} instances: Rust {rust_secs:.1f}s total; W_ext mismatches {bad}; worst len/OPT {worst:.3f}; rounds {dict(sorted(rounds.items()))}")
