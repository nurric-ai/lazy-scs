"""Search for instances that hit rare branches, using the fast Rust binary."""
import random, subprocess, sys, os, json
BIN = os.path.join(os.environ["TEMP"], "scsrs", "release", "scs-rs.exe")
RARE = ["requests", "host_requests", "individual_link", "below_baseline", "self_cycle",
        "upper_round_trip", "short_contact_link", "opened_cycle", "band_join"]

def gen(rng):
    kind = rng.randrange(4)
    out = []
    for _ in range(rng.randint(2, 9)):
        if kind == 0:      # nested periodic texts with flanks
            small = rng.choice(["a", "ab", "aab", "abb", "abc", "aabb", "abab"])
            prim = small * rng.randint(1, 6) + rng.choice(["", "c", "b", "cc", "bc", "d", "cd", "ac"])
        elif kind == 1:    # random primitive words over a tiny alphabet
            prim = "".join(rng.choice("ab" if rng.random() < .7 else "abc") for _ in range(rng.randint(1, 9)))
        elif kind == 2:    # mixtures of two periods
            u, v = rng.choice(["a", "ab", "aab"]), rng.choice(["b", "ba", "abb", "bba"])
            prim = u * rng.randint(1, 5) + v * rng.randint(1, 5)
        else:              # squares and near-squares
            w = "".join(rng.choice("ab") for _ in range(rng.randint(2, 7)))
            prim = w * rng.randint(2, 4) + rng.choice(["", "a", "b", "c"])
        L = rng.randint(2, 40)
        off = rng.randrange(len(prim))
        s = (prim * (L // len(prim) + 3))[off: off + L]
        if rng.random() < 0.15:
            i = rng.randrange(len(s)); s = s[:i] + rng.choice("abc") + s[i + 1:]
        out.append(s)
    return out

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
N = int(sys.argv[1]); chunk = 2000
TRIG = sys.argv[3] if len(sys.argv) > 3 else "triggers.json"
found, hits = {}, {k: 0 for k in RARE}
for c0 in range(0, N, chunk):
    insts = [gen(rng) for _ in range(min(chunk, N - c0))]
    open("hunt_batch.txt", "w").write("\n---\n".join("\n".join(i) for i in insts) + "\n")
    out = subprocess.run([BIN, "--batch", "hunt_batch.txt"], capture_output=True, text=True)
    if out.returncode:
        print("CRASH", out.stderr[-1500:]); break
    lines = out.stdout.splitlines()
    for i, inst in enumerate(insts):
        head = lines[2 * i].split()
        stats = dict(kv.split("=") for kv in head[4].split(",")) if len(head) > 4 else {}
        for k in RARE:
            if k in stats:
                hits[k] += 1
                size = sum(map(len, inst))
                if k not in found or size < sum(map(len, found[k])):
                    found[k] = inst
    print(c0 + len(insts), {k: v for k, v in hits.items() if v}, flush=True)
    json.dump(found, open(TRIG, "w"), indent=1)
print(N, "instances; hits:", hits)
json.dump(found, open(TRIG, "w"), indent=1)
for k, v in found.items():
    print(k, sum(map(len, v)), v)
