"""Baselines: greedy merge, exact subset DP, and an output trimmer."""
from scs2 import reduce_inputs


def overlap(a, b):
    """Longest proper suffix of a that is a prefix of b (prefix function)."""
    s = b + "\0" + a
    pi = [0] * len(s)
    for i in range(1, len(s)):
        k = pi[i - 1]
        while k and s[i] != s[k]:
            k = pi[k - 1]
        if s[i] == s[k]:
            k += 1
        pi[i] = k
    return min(pi[-1], len(a) - 1, len(b) - 1)


def greedy(strings):
    """Classic greedy: repeatedly merge the pair with the largest overlap."""
    S = reduce_inputs(strings)
    if not S:
        return ""
    alive = list(range(len(S)))
    ov = {(i, j): overlap(S[i], S[j]) for i in alive for j in alive if i != j}
    while len(alive) > 1:
        i, j = max(((i, j) for i in alive for j in alive if i != j),
                   key=lambda ij: (ov[ij], -ij[0], -ij[1]))
        S.append(S[i] + S[j][ov[i, j]:])
        n = len(S) - 1
        alive = [x for x in alive if x not in (i, j)]
        for x in alive:
            ov[n, x] = overlap(S[n], S[x])
            ov[x, n] = overlap(S[x], S[n])
        alive.append(n)
        # a merged string may now contain another alive string
        alive = [x for x in alive if x == n or S[x] not in S[n]]
    return S[alive[0]]


def exact(strings, max_n=13):
    """Optimal superstring by DP over subsets (Held-Karp on overlaps)."""
    S = reduce_inputs(strings)
    n = len(S)
    if n == 0:
        return ""
    assert n <= max_n
    ov = [[overlap(S[i], S[j]) if i != j else 0 for j in range(n)] for i in range(n)]
    INF = float("inf")
    dp = [[INF] * n for _ in range(1 << n)]
    par = [[-1] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = len(S[i])
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if cur == INF:
                continue
            for nx in range(n):
                if mask >> nx & 1:
                    continue
                v = cur + len(S[nx]) - ov[last][nx]
                nm = mask | 1 << nx
                if v < dp[nm][nx]:
                    dp[nm][nx] = v
                    par[nm][nx] = last
    full = (1 << n) - 1
    last = min(range(n), key=lambda i: dp[full][i])
    order, mask = [], full
    while last != -1:
        order.append(last)
        pl = par[mask][last]
        mask ^= 1 << last
        last = pl
    order.reverse()
    out = S[order[0]]
    for a, b in zip(order, order[1:]):
        out += S[b][ov[a][b]:]
    return out


def trim(T, strings):
    """Drop leading/trailing letters while every input still occurs (never lengthens)."""
    req = [s for s in strings if s]
    while T and all(s in T[1:] for s in req):
        T = T[1:]
    while T and all(s in T[:-1] for s in req):
        T = T[:-1]
    return T
