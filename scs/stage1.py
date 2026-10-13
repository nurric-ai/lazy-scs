"""Stage 1 go/no-go: which periodic rules (2.6) actually matter?

Rule (v, p, k): v has period p < |v|, w = v minus k periods.
  square rule : |v| >= 2p   (v is a repetition -> lives inside a run; few of them)
  border rule : |v| <  2p   (v merely has a border; present almost everywhere)

For each word w we compare
  base   = extension bounds (2.5) and the required-string bound
  sq     = best bound from fired square rules
  border = best bound from fired border rules
If border <= max(base, sq) at every w, dropping all border rules leaves every
count m(s) unchanged (induction on decreasing length), so W and the whole
algorithm are identical. We also log the rules the connection phase relies on
(Section 5.3) and re-run the full pipeline with square rules only.

    python stage1.py
"""
import random
import sys
from collections import Counter

from scs2 import Counts, reduce_inputs, two_approx
import fuzz
import bench
import shibata


def is_square(v, p, k):
    return len(v) >= 2 * p


def sources(rng):
    for _ in range(3000):
        yield "fuzz-mixed", fuzz.gen(rng)
    for _ in range(400):
        yield "random-4", bench.fam_random(rng, rng.randint(3, 12), 12, "abcd")
    for _ in range(400):
        yield "reads", bench.fam_reads(rng, rng.randint(5, 30), 20, 120)
    for _ in range(400):
        yield "periodic", bench.fam_periodic(rng, rng.randint(3, 8), 16)
    for _ in range(300):
        out = []
        for _ in range(rng.randint(2, 6)):
            small = rng.choice(["a", "ab", "aab", "abc", "abb", "aabb"])
            prim = small * rng.randint(1, 4) + rng.choice(["", "c", "b", "bc", "d"])
            L = rng.randint(2, 24)
            off = rng.randrange(len(prim))
            out.append((prim * (L // len(prim) + 3))[off: off + L])
        yield "nested", out
    for k, t in [(10, 2), (10, 4), (12, 3)]:
        yield "shibata", shibata.instance(k, t)[0]


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    agg = Counter()
    per_src = Counter()
    examples = []
    used_exp = Counter()
    for src, inst in sources(rng):
        req = reduce_inputs(inst)
        if not req:
            continue
        events = Counter()

        def log(w, base, fired):
            sq = max([b for v, p, k, b in fired if is_square(v, p, k)], default=0)
            bo = max([b for v, p, k, b in fired if not is_square(v, p, k)], default=0)
            events["words"] += 1
            if max(sq, bo) > base:
                events["rules_matter"] += 1
            if sq > base and sq >= bo:
                events["square_binds"] += 1
            if bo > max(base, sq):
                events["border_binds"] += 1
                if len(examples) < 5:
                    examples.append((src, inst, w, base, sq, bo,
                                     [(v, p, k, b) for v, p, k, b in fired if not is_square(v, p, k)][:3]))

        C = Counts(req, log=log)
        Csq = Counts(req, rule_filter=is_square)
        agg["instances"] += 1
        agg["words"] += events["words"]
        agg["rules_total"] += C.n_rules
        agg["rules_square"] += Csq.n_rules
        for key in ("rules_matter", "square_binds", "border_binds"):
            if events[key]:
                agg["inst_" + key] += 1
                per_src[src, key] += 1
        per_src[src, "n"] += 1
        if C.m != Csq.m:
            agg["inst_m_differs"] += 1

        # full pipeline with square rules only must still work
        try:
            r = two_approx(inst, rule_filter=is_square)
            assert all(s in r.superstring for s in inst) and len(r.superstring) <= 2 * r.W
            agg["sq_pipeline_ok"] += 1
        except AssertionError as e:
            agg["sq_pipeline_FAIL"] += 1
            if len(examples) < 8:
                examples.append(("PIPELINE FAIL", inst, repr(e)))
        rf = two_approx(inst)
        for v, p, k in rf.used_rules:
            used_exp["square" if len(v) >= 2 * p else "border"] += 1

    print("== totals ==")
    for k, v in sorted(agg.items()):
        print(f"  {k}: {v}")
    print(f"  rule count: square rules are {100 * agg['rules_square'] / agg['rules_total']:.1f}% of all rules")
    print("== rules relied on by the connection phase (5.3) ==", dict(used_exp))
    print("== per source: instances where [rules matter / square binds / border binds] ==")
    for src in sorted({s for s, _ in per_src}):
        n = per_src[src, "n"]
        print(f"  {src:12s} n={n:5d}  matter={per_src[src, 'rules_matter']:5d}  "
              f"square={per_src[src, 'square_binds']:5d}  border={per_src[src, 'border_binds']:5d}")
    for ex in examples:
        print("EXAMPLE:", ex)


if __name__ == "__main__":
    main()
