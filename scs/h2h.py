"""Head-to-head: our Rust implementation vs the JavaScript (bkarak/shortest-superstring-js)
and C++ (spele1100/scs-2approx) implementations of openai/math #128.

Inputs: error-free 150 bp reads at 10x coverage from prefixes of the E. coli K-12
genome (NC_000913.3). Each run gets a time limit; once an implementation times out
or runs out of memory it is skipped for larger sizes.

    python h2h.py <ecoli.fa> <js index.js> <cpp scs.exe> <rust scs-rs.exe> [limit_s]
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import time

import psutil


def run(cmd, limit):
    t0 = time.time()
    p = psutil.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    peak = 0
    while p.poll() is None:
        try:
            peak = max(peak, p.memory_info().rss)
        except psutil.Error:
            pass
        if time.time() - t0 > limit:
            p.kill()
            return None, time.time() - t0, peak, "timeout"
        time.sleep(0.05)
    out, err = p.communicate()
    if p.returncode != 0:
        return None, time.time() - t0, peak, f"exit {p.returncode}: {err.strip()[-120:]}"
    return out, time.time() - t0, peak, "ok"


def main():
    fa, js, cpp, rs = sys.argv[1:5]
    limit = float(sys.argv[5]) if len(sys.argv) > 5 else 300
    genome = "".join(l.strip() for l in open(fa) if not l.startswith(">"))
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "h2h_run_js.js")
    only = set(sys.argv[6].split(",")) if len(sys.argv) > 6 else {"js", "cpp", "rust"}
    alive = {k: k in only for k in ("js", "cpp", "rust")}
    rows = []
    tmp = tempfile.mkdtemp()
    for G in [1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000, len(genome)]:
        g = genome[:G]
        rng = random.Random(G)
        reads = [g[i:i + 150] for i in (rng.randrange(G - 149) for _ in range(G * 10 // 150))]
        f1 = os.path.join(tmp, "reads.txt")
        f2 = os.path.join(tmp, "reads_cpp.txt")
        open(f1, "w", newline="\n").write("\n".join(reads) + "\n")
        open(f2, "w", newline="\n").write(f"{len(reads)}\n" + "\n".join(reads) + "\n")
        row = {"genome": G, "reads": len(reads), "chars": 150 * len(reads)}
        for name, cmd, parse in [
            ("js", ["node", "--max-old-space-size=16000", runner, f1, js], lambda o: (json.loads(o)["length"], json.loads(o)["lowerBound"])),
            ("cpp", [cpp, f2, "--json"], lambda o: (json.loads(o)["len"], json.loads(o)["W"])),
            ("rust", [rs, f1], lambda o: (json.loads(o)["length"], json.loads(o)["W"])),
        ]:
            if not alive[name]:
                row[name] = "skipped" if name in only else "not run"
                continue
            out, secs, peak, status = run(cmd, limit)
            if status != "ok":
                alive[name] = False
                row[name] = f"{status} after {secs:.0f}s, {peak / 2**20:.0f} MB"
            else:
                ln, w = parse(out)
                row[name] = f"{secs:.2f}s, {peak / 2**20:.0f} MB, len {ln}, W {w}"
        rows.append(row)
        print(json.dumps(row), flush=True)
    json.dump(rows, open(f"h2h_results_{'_'.join(sorted(only))}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
