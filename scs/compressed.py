"""Stage 2a: extension-only counts m(s) without enumerating substrings.

Generalized suffix automaton over the (reduced) inputs. A state X holds the
substrings with lengths (len(link X), len X]; they share right extensions
(transitions), and every member except the longest has exactly one left
extension (the next longer member). So with only the extension rules (2.5)
and the required-string rule, m restricted to X is a nonincreasing step
function of length:

  m_X(hi) = max(sum_{Z: link Z = X} m_Z(minlen Z),  R_X(hi),  [X is an input])
  m_X(l)  = max(m_X(l+1), R_X(l))                         for l < hi
  R_X(l)  = sum_c m_{delta(X,c)}(l+1)

States are processed by decreasing len, which respects every dependency.
This is the count table that the lazy-rule algorithm starts from.
"""
from bisect import bisect_right

from scs2 import reduce_inputs


class SAM:
    def __init__(self, strings):
        self.next, self.link, self.len, self.pos = [{}], [-1], [0], [(0, 0)]
        self.strings = strings
        for si, s in enumerate(strings):
            last = 0
            for i, c in enumerate(s):
                self._cur_pos = (si, i + 1)
                last = self._extend(last, c)

    def _new(self, ln, nxt=None, link=-1, pos=None):
        self.pos.append(pos if pos is not None else self._cur_pos)
        self.next.append(dict(nxt) if nxt else {})
        self.link.append(link)
        self.len.append(ln)
        return len(self.len) - 1

    def _split(self, p, c, q):
        clone = self._new(self.len[p] + 1, self.next[q], self.link[q], self.pos[q])
        while p != -1 and self.next[p].get(c) == q:
            self.next[p][c] = clone
            p = self.link[p]
        self.link[q] = clone
        return clone

    def _extend(self, last, c):
        nx, ln = self.next, self.len
        if c in nx[last]:                       # generalized SAM: path already exists
            q = nx[last][c]
            return q if ln[last] + 1 == ln[q] else self._split(last, c, q)
        cur = self._new(ln[last] + 1)
        p = last
        while p != -1 and c not in nx[p]:
            nx[p][c] = cur
            p = self.link[p]
        if p == -1:
            self.link[cur] = 0
        else:
            q = nx[p][c]
            self.link[cur] = q if ln[p] + 1 == ln[q] else self._split(p, c, q)
        return cur

    def word(self, X, l):
        si, e = self.pos[X]
        return self.strings[si][e - l: e]

    def __contains__(self, s):
        v = 0
        for c in s:
            v = self.next[v].get(c)
            if v is None:
                return False
        return True

    def state_of(self, s):
        v = 0
        for c in s:
            v = self.next[v][c]
        return v


def ev(f, l):
    """Evaluate a step function [(start, value), ...] at length l."""
    return f[bisect_right(f, (l, float("inf"))) - 1][1]


class CompressedCounts:
    def __init__(self, strings, texts=(), max_passes=1000):
        """texts: enabled periodic texts (p, least rotation) whose rules (2.6) apply.
        Rules are applied as point lower bounds ("bumps") and iterated to a fixpoint;
        each bump is valid for any count vector bounded by the true occurrence counts."""
        req = reduce_inputs(strings)
        self.sam = SAM(req)
        self.req = req
        self.bumps = {}
        self._compute()
        if not texts:
            return
        rules = self._rules(texts)
        self.passes = 0
        for self.passes in range(1, max_passes + 1):
            changed = False
            cache = {}
            for vX, vl, wX, wl, P, k in rules:
                if ev(self.f[vX], vl) > self._B(P, vl, cache):
                    need = self._B(P, wl, cache) + k + 1
                    if ev(self.f[wX], wl) < need:
                        self.bumps[wX, wl] = need
                        changed = True
            if not changed:
                break
            self._compute()

    def _rules(self, texts):
        """All (v, p, k) with v in V on an enabled periodic text, as (state, len) pairs."""
        S, out = self.sam, []
        for p, G in texts:
            for r in range(p):
                P = G[r:] + G[:r]
                path, v, i = [0], 0, 0
                while True:
                    v = S.next[v].get(P[i % p])
                    if v is None or i > 4 * max(map(len, self.req)):
                        break
                    path.append(v)
                    i += 1
                for n in range(p + 1, len(path)):
                    for k in range(1, (n - 1) // p + 1):
                        out.append((path[n], n, path[n - k * p], n - k * p, P, k))
        return out

    def _m_word(self, s):
        v = 0
        for c in s:
            v = self.sam.next[v].get(c)
            if v is None:
                return 0
        return ev(self.f[v], len(s)) if s else 0

    def _B(self, P, slen, cache):
        """Blocking functional (2.3) via the automaton (same enumeration as the reference)."""
        key = (P, slen)
        if key in cache:
            return cache[key]
        S, p, Lm = self.sam, len(P), max(map(len, self.req))
        off = Lm + 1
        span = off + slen + Lm + 2
        st = (-off) % p
        Astr = (P * ((st + span) // p + 2))[st: st + span]
        sigma = sorted({c for s in self.req for c in s})
        total, a = 0, 0
        while a + slen <= Lm and Astr[off - a: off + slen] in S:
            c = 0
            while a + slen + c <= Lm:
                M = Astr[off - a: off + slen + c]
                if M not in S:
                    break
                xbad, ybad = Astr[off - a - 1], Astr[off + slen + c]
                for x in sigma:
                    if x == xbad or (x + M) not in S:
                        continue
                    for y in sigma:
                        if y != ybad:
                            total += self._m_word(x + M + y)
                c += 1
            a += 1
        cache[key] = total
        return total

    def _compute(self):
        S = self.sam
        req = self.req
        n = len(S.len)
        self.is_req = [False] * n
        for s in req:
            self.is_req[S.state_of(s)] = True
        children_sum = [0] * n
        bump_at = {}
        for (Y, l) in self.bumps:
            bump_at.setdefault(Y, []).append(l)
        self.f = f = [None] * n     # f[X] = step function on [lo, hi]
        self.R = R = [None] * n
        for X in sorted(range(1, n), key=lambda x: -S.len[x]):
            lo, hi = S.len[S.link[X]] + 1, S.len[X]
            # R_X(l) = sum_c m_Y(l + 1): breakpoints are Y's starts shifted by one
            ys = list(S.next[X].values())
            pts = {lo}
            for Y in ys:
                for st, _ in f[Y]:
                    if lo <= st - 1 <= hi:
                        pts.add(st - 1)
            pts = sorted(pts)
            r = [(l, sum(ev(f[Y], l + 1) for Y in ys)) for l in pts]
            R[X] = r
            # running max from hi downward; bumps are point lower bounds at (X, l)
            cur = max(children_sum[X], 1 if self.is_req[X] else 0)
            bl = sorted(bump_at.get(X, ()))
            if bl:
                cut = sorted(set(pts) | set(bl) | {l + 1 for l in bl if l < hi})
                r = [(l, ev(r, l)) for l in cut]
            segs = []
            for l, val in reversed(r):
                val = max(val, self.bumps.get((X, l), 0)) if bl else val
                cur = max(cur, val)
                segs.append((l, cur))
            segs.reverse()
            merged = []
            for l, val in segs:
                if not merged or merged[-1][1] != val:
                    merged.append((l, val))
            f[X] = merged
            children_sum[S.link[X]] += merged[0][1]   # m_X(minlen X) feeds its link's left sum
        self.children_sum = children_sum
        self.W = sum(ev(f[Y], 1) for Y in S.next[0].values())
        self.breakpoints = sum(len(x) for x in f[1:])
        self.states = n - 1

    def m(self, s):
        return ev(self.f[self.sam.state_of(s)], len(s))

    def edges(self):
        """Nonzero up/down multiplicities keyed by word, enumerated segment-wise:
        total work O(states + W), never touching zero-edge substrings."""
        from collections import Counter
        S, up, down = self.sam, Counter(), Counter()
        for X in range(1, len(S.len)):
            lo, hi, fx, rx = S.len[S.link[X]] + 1, S.len[X], self.f[X], self.R[X]
            u_hi = ev(fx, hi) - self.children_sum[X]
            if u_hi:
                up[S.word(X, hi)] = u_hi
            for st, val in fx[1:]:               # f drops at st: u(st - 1) = f(st - 1) - f(st)
                up[S.word(X, st - 1)] = ev(fx, st - 1) - val
            cuts = sorted({l for l, _ in fx} | {l for l, _ in rx})
            for i, a in enumerate(cuts):
                b = cuts[i + 1] - 1 if i + 1 < len(cuts) else hi
                d = ev(fx, a) - ev(rx, a)
                if d:
                    for l in range(a, b + 1):
                        down[S.word(X, l)] = d
        return up, down
