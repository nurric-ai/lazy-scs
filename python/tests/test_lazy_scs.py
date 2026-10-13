"""Tests for the installed lazy_scs package: pytest python/tests"""
import itertools
import random
import subprocess
import sys

import pytest

import lazy_scs


def brute_force_opt(strings):
    """Exact SCS length for tiny inputs (permutations with maximal overlaps)."""
    strings = [s for s in strings if not any(s != t and s in t for t in strings)]
    strings = list(dict.fromkeys(strings))
    best = None
    for perm in itertools.permutations(strings):
        cur = perm[0]
        for s in perm[1:]:
            k = next(k for k in range(min(len(cur), len(s)), -1, -1) if k == 0 or cur.endswith(s[:k]))
            cur += s[k:]
        best = len(cur) if best is None else min(best, len(cur))
    return best or 0


def check(strings, r):
    for s in strings:
        assert s in r.superstring
    assert r.length <= 2 * r.lower_bound


def test_example():
    r = lazy_scs.solve(["GATTACA", "TACAGG", "AGGCT"])
    check(["GATTACA", "TACAGG", "AGGCT"], r)
    assert isinstance(r.superstring, str)
    assert (r.superstring, r.lower_bound, r.ratio) == ("GATTACAGGCT", 11, 1.0)


def test_bytes_in_bytes_out():
    r = lazy_scs.solve([b"abc", b"bcd", b"cde"])
    assert isinstance(r.superstring, bytes)
    check([b"abc", b"bcd", b"cde"], r)


def test_edge_cases():
    assert lazy_scs.shortest_superstring([]) == ""
    assert lazy_scs.shortest_superstring(["abc"]) == "abc"
    assert lazy_scs.shortest_superstring(["abc", "abc", "b"]) == "abc"
    r = lazy_scs.solve(["aaaa", "aaa", "aa"])           # periodic input
    check(["aaaa", "aaa", "aa"], r)


def test_mixed_types_rejected():
    with pytest.raises(TypeError):
        lazy_scs.solve(["abc", b"abc"])


@pytest.mark.parametrize("seed", range(300))
def test_random_against_exact_optimum(seed):
    rng = random.Random(seed)
    alphabet = "ab" if seed % 3 == 0 else "abc"
    strings = ["".join(rng.choice(alphabet) for _ in range(rng.randint(1, 6))) for _ in range(rng.randint(1, 6))]
    r = lazy_scs.solve(strings, verify=True)
    check(strings, r)
    opt = brute_force_opt(strings)
    assert r.lower_bound <= opt <= r.length <= 2 * opt


def test_reads_from_a_genome():
    rng = random.Random(7)
    genome = "".join(rng.choice("ACGT") for _ in range(20000))
    reads = [genome[i:i + 150] for i in sorted(rng.sample(range(len(genome) - 150), 1500))]
    r = lazy_scs.solve(reads, verify=True)
    check(reads, r)
    assert r.ratio < 1.01


def test_cli(tmp_path):
    p = tmp_path / "in.txt"
    p.write_text(">r1\nACGTT\n>r2\nGTTAC\n>r3\nTACCA\n")
    out = tmp_path / "out.txt"
    subprocess.run([sys.executable, "-m", "lazy_scs", str(p), "-o", str(out), "--verify"], check=True)
    sup = out.read_text()
    assert all(s in sup for s in ["ACGTT", "GTTAC", "TACCA"])
