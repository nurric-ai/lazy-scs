"""Command line: lazy-scs INPUT [-o OUTPUT] [--verify] [--json]

INPUT has one string per line, or FASTA records (any line starting with '>'); '-' reads stdin.
Prints the superstring (or writes it to OUTPUT) and a summary line on stderr.
"""
from __future__ import annotations

import argparse
import json
import sys
import time

from . import __version__, solve


def read_strings(text: bytes) -> list[bytes]:
    lines = text.splitlines()
    if any(l.startswith(b">") for l in lines):
        out, cur = [], []
        for l in lines:
            if l.startswith(b">"):
                if cur:
                    out.append(b"".join(cur))
                cur = []
            else:
                cur.append(l.strip())
        if cur:
            out.append(b"".join(cur))
        return [s for s in out if s]
    return [l.rstrip(b"\r\n") for l in lines if l.rstrip(b"\r\n")]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="lazy-scs", description="Factor-2 approximate shortest common superstring.")
    ap.add_argument("input", help="file with one string per line or FASTA; '-' for stdin")
    ap.add_argument("-o", "--output", help="write the superstring here instead of stdout")
    ap.add_argument("--verify", action="store_true", help="independently re-check the lower bound")
    ap.add_argument("--json", action="store_true", help="print the summary as JSON on stdout (no superstring)")
    ap.add_argument("--version", action="version", version=f"lazy-scs {__version__}")
    a = ap.parse_args(argv)

    raw = sys.stdin.buffer.read() if a.input == "-" else open(a.input, "rb").read()
    strings = read_strings(raw)
    t0 = time.perf_counter()
    r = solve(strings, verify=a.verify)
    summary = {"strings": len(strings), "input_chars": sum(map(len, strings)), "seconds": round(time.perf_counter() - t0, 3),
               "length": r.length, "lower_bound": r.lower_bound, "ratio": round(r.ratio, 6), "rounds": r.rounds}
    if a.verify:
        summary["verified"] = True
    if a.output:
        with open(a.output, "wb") as f:
            f.write(r.superstring)
    if a.json:
        print(json.dumps(summary))
    else:
        if not a.output:
            sys.stdout.buffer.write(r.superstring + b"\n")
        print(json.dumps(summary), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
