"""Greedy beats 2x OPT: Shibata's counterexample (arXiv:2609.01365) vs the 2-approximation.

Builds S_{k,t}, replays the paper's greedy run (checking that every merge is a
legal greedy step, i.e. has maximum overlap at that moment), and compares it with
the factor-2 algorithm on the same instance.

    python shibata.py
"""
import time

from scs2 import two_approx
from baselines import overlap, greedy, trim


def X(n):
    return "yx" * (n // 2) if n % 2 == 0 else "x" + "yx" * ((n - 1) // 2)


def csub(T, j, l):            # CSubstr_l(T, j), 1-indexed
    rep = T * (l // len(T) + 3)
    return rep[j - 1: j - 1 + l]


def instance(k, t):
    s = k // 2
    U, V, W, S = {}, {}, {}, []
    for i in range(1, t + 1):
        b = chr(0x3B1 + i - 1)            # b_i = distinct Greek letters
        U[i] = "xa" + X(s - 3) + b + "xa" + X(s - 3) + "y"
        V[i] = "xa" + X(s - 3) + b
        W[i] = "xa" + X(s - 1) + b
        assert len(U[i]) == k and len(V[i]) == s and len(W[i]) == s + 2
        S += [csub(U[i], j, k) for j in range(1, k + 1)]
        S += [csub(V[i], j, k) for j in range(1, s + 1)]
        S += [csub(W[i], j, k) for j in range(1, s + 3)]
    assert len(set(S)) == len(S)
    return S, U, V, W


def scripted_greedy(S, U, V, W, k, t):
    """Replay Shibata's merge order; every merge must have the current max overlap."""
    s = k // 2
    pool = set(S)
    ov = {(a, b): overlap(a, b) for a in pool for b in pool if a != b}

    def merge(a, b):
        m = max(ov.values())
        o = overlap(a, b)
        assert o == m, f"illegal greedy step: overlap {o} < max {m}"
        z = a + b[o:]
        pool.discard(a); pool.discard(b)
        for key in [key for key in ov if a in key or b in key]:
            del ov[key]
        for c in pool:
            ov[z, c] = overlap(z, c)
            ov[c, z] = overlap(c, z)
        pool.add(z)
        return z

    def chain(T, start, count):
        cur = csub(T, start, k)
        for r in range(1, count):
            cur = merge(cur, csub(T, (start - 1 + r) % len(T) + 1, k))
        return cur

    hatU = {i: chain(U[i], 4, k) for i in range(1, t + 1)}
    hatV = {i: chain(V[i], 2, s) for i in range(1, t + 1)}
    hatW = {i: chain(W[i], 6, s + 2) for i in range(1, t + 1)}
    M = {i: merge(hatV[i], hatU[i]) for i in range(1, t + 1)}
    cur = M[1]
    for i in range(2, t + 1):
        cur = merge(cur, M[i])
    for i in range(1, t + 1):
        cur = merge(cur, hatW[i])
    assert len(pool) == 1
    return cur


if __name__ == "__main__":
    print("| k | t | strings | OPT | greedy (bad ties) | ratio | greedy (our ties) | 2apx | 2apx+trim | ratio | 2apx time |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for k, t in [(10, 5), (10, 20), (14, 10), (20, 10)]:
        S, U, V, W = instance(k, t)
        s = k // 2
        opt = t * (4 * s + 2) + k - 1                 # Theorem 10
        bad = scripted_greedy(S, U, V, W, k, t)
        assert len(bad) == t * (9 * s + 1) + 1 and all(x in bad for x in S)   # Theorem 13
        g = greedy(S)
        t0 = time.time()
        r = two_approx(S)
        dt = time.time() - t0
        tr = trim(r.superstring, S)
        assert r.W <= opt and all(x in tr for x in S)
        print(f"| {k} | {t} | {len(S)} | {opt} | {len(bad)} | **{len(bad)/opt:.3f}** | {len(g)} | "
              f"{len(r.superstring)} | {len(tr)} | **{len(tr)/opt:.3f}** | {dt:.1f}s |")
