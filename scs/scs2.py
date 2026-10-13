"""
Polynomial-time 2-approximation for Shortest Common Superstring.

Implements openai/math result 128, "A Polynomial-Time 2-Approximation for
Shortest Common Superstring" (OpenAI, Sept 24 2026). Section/lemma numbers in
comments refer to that paper.

    from scs2 import two_approx
    res = two_approx(["cde", "abc", "efg"])
    res.superstring, res.W          # 'abcdefg', 7   (len <= 2W <= 2 OPT)

Reference implementation: favours following the paper over speed. Time and
memory grow with the number of distinct substrings, i.e. O(L^2) words.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field

CHECK = True  # assert the paper's invariants while running


def reduce_inputs(strings):
    """Drop empties, duplicates and strings contained in another (OPT-preserving)."""
    by_len = sorted(set(s for s in strings if s), key=len, reverse=True)
    kept = []
    for s in by_len:
        if not any(s in k for k in kept):
            kept.append(s)
    return sorted(kept)


# --------------------------------------------------------------------------
# Section 2: forced occurrence counts m(s) and the balanced base graph
# --------------------------------------------------------------------------
class Counts:  # reference: one entry per distinct substring
    def __init__(self, req, rule_filter=None, log=None):
        """rule_filter(v, p, k) -> bool keeps a subset of rules (2.6); any subset still
        gives a valid lower bound by Lemma 2.1. log(w, base, [(v,p,k,bound)]) observes."""
        V = {""}
        for s in req:
            for i in range(len(s)):
                for j in range(i + 1, len(s) + 1):
                    V.add(s[i:j])
        self.V = V
        self.req = set(req)
        lext, rext, brackets = defaultdict(list), defaultdict(list), defaultdict(list)
        for r in V:
            if r:
                lext[r[1:]].append(r)   # r = c + s
                rext[r[:-1]].append(r)  # r = s + c
            if len(r) >= 2:
                brackets[r[1:-1]].append(r)
        self.lext, self.rext, self.brackets = lext, rext, brackets

        # Rules (2.6): v in V with period p < |v|; w = v minus k periods.
        rules = defaultdict(list)
        for v in V:
            n = len(v)
            for p in range(1, n):
                if v[p:] == v[:-p]:
                    for k in range(1, (n - 1) // p + 1):
                        if rule_filter is None or rule_filter(v, p, k):
                            rules[v[: n - k * p]].append((v, p, k))
        self.n_rules = sum(map(len, rules.values()))

        self.m = m = {}
        self._bcache = {}
        self.maxlen = max(map(len, req))
        for s in sorted(V - {""}, key=len, reverse=True):
            val = max(sum(m.get(r, 0) for r in lext[s]), sum(m.get(r, 0) for r in rext[s]))
            if s in self.req:
                val = max(val, 1)
            fired = []
            for v, p, k in rules.get(s, ()):
                P = v[:p]
                if m[v] > self.B(P, len(v)):
                    fired.append((v, p, k, self.B(P, len(s)) + k + 1))
            if log is not None:
                log(s, val, fired)
            m[s] = max([val] + [f[3] for f in fired])
        self.W = sum(m.get(c, 0) for c in {s[0] for s in V if s})

    def B(self, P, slen):
        """Blocking functional (2.3) for s = A[0:slen], A = P repeated, aligned at 0."""
        p = next(d for d in range(1, len(P) + 1) if len(P) % d == 0 and P == P[:d] * (len(P) // d))
        P = P[:p]                                   # same text, canonical key
        key = (P, slen)
        if key in self._bcache:
            return self._bcache[key]
        V, m, Lm = self.V, self.m, self.maxlen
        off = Lm + 1                                # Astr[off + i] = A(i)
        span = off + slen + Lm + 2
        st = (-off) % p
        Astr = (P * ((st + span) // p + 2))[st: st + span]
        total, a = 0, 0
        while a + slen <= Lm and Astr[off - a: off + slen] in V:
            c = 0
            while a + slen + c <= Lm:
                M = Astr[off - a: off + slen + c]
                if M not in V:
                    break
                xbad, ybad = Astr[off - a - 1], Astr[off + slen + c]
                for r in self.brackets.get(M, ()):
                    if r[0] != xbad and r[-1] != ybad:
                        total += m.get(r, 0)
                c += 1
            a += 1
        self._bcache[key] = total
        return total

    def edges(self):
        """(2.8): up-edge multiplicity u(s) into s, down-edge multiplicity d(s) out of s."""
        m, up, down = self.m, Counter(), Counter()
        for s, ms in m.items():
            u = ms - sum(m.get(r, 0) for r in self.lext[s])
            d = ms - sum(m.get(r, 0) for r in self.rext[s])
            if CHECK:
                assert u >= 0 and d >= 0
            if u:
                up[s] = u
            if d:
                down[s] = d
        return up, down


# --------------------------------------------------------------------------
# Lemma 4.1: following ordered windows; walks are recorded as edge multisets
# --------------------------------------------------------------------------
class Walks:
    def __init__(self, V):
        self.up, self.down, self.cost, self.V = Counter(), Counter(), 0, V

    def follow(self, text, windows):
        word = lambda a, b: "".join(text(i) for i in range(a, b))
        if CHECK:
            for a, b in windows:
                assert a <= b and word(a, b) in self.V, ("window not a vertex", a, b)
            assert word(*windows[0]) == word(*windows[-1]), "walk not closed"
        for (a, b), (c, d) in zip(windows, windows[1:]):
            assert a <= c and b <= d, ("unordered windows", (a, b), (c, d))
            if c <= b:
                for x in range(a, c):
                    self.down[word(x, b)] += 1
            else:
                for x in range(a, b):
                    self.down[word(x, b)] += 1
                for i in range(b, c):
                    ch = text(i)
                    self.up[ch] += 1
                    self.down[ch] += 1
            for e in range(max(b, c) + 1, d + 1):
                self.up[word(c, e)] += 1
            self.cost += d - b

    def round_trip(self, s):
        """epsilon -> s -> epsilon, cost |s|."""
        for i in range(1, len(s) + 1):
            self.up[s[:i]] += 1
            self.down[s[i - 1:]] += 1
        self.cost += len(s)


# --------------------------------------------------------------------------
# Section 3: closed walks -> ordered periodic layers
# --------------------------------------------------------------------------
def closed_walks(up, down):
    """Decompose the balanced base graph into closed walks (lists of steps)."""
    out = defaultdict(list)  # vertex -> [ [kind, word, remaining] ]
    for s in sorted(up):
        out[s[:-1]].append(["u", s, up[s]])
    for s in sorted(down):
        out[s].append(["d", s, down[s]])
    ptr = defaultdict(int)

    def take(v):
        lst = out[v]
        while ptr[v] < len(lst) and lst[ptr[v]][2] == 0:
            ptr[v] += 1
        if ptr[v] == len(lst):
            return None
        e = lst[ptr[v]]
        e[2] -= 1
        return e

    walks = []
    for start in sorted(out, key=lambda v: (len(v), v)):
        while True:
            e = take(start)
            if e is None:
                break
            steps, v = [], start
            while True:
                kind, s, _ = e
                if kind == "u":
                    steps.append(("u", s[-1]))
                    v = s
                else:
                    steps.append(("d", None))
                    v = s[1:]
                if v == start:
                    break
                e = take(v)
                assert e is not None, "graph not balanced"
            walks.append((start, steps))
    return walks


def least_rotation(s):
    return min(s[i:] + s[:i] for i in range(len(s)))


@dataclass
class Group:
    p: int
    G: str              # least rotation of the primitive text; A(i) = G[i % p]
    t: int = 0          # distinguished position: start of the greatest rotation
    layers: list = field(default_factory=list)  # global layer ids, ascending

    def A(self, i):
        return self.G[i % self.p]

    def word(self, a, b):
        return "".join(self.G[i % self.p] for i in range(a, b))


@dataclass
class Layer:
    gid: int
    zb: list            # z(x) for x in [0, p)

    def z(self, x, p):
        return self.zb[x % p] + p * (x // p)


def build_layers(walks):
    groups_by_key, groups, layers = {}, [], []
    gwalks = defaultdict(list)  # gid -> list of (P, Zg) with Zg over [0,P) in group coords
    for start, steps in walks:
        U = "".join(ch for k, ch in steps if k == "u")
        P = len(U)
        p = next(q for q in range(1, P + 1) if P % q == 0 and U == U[q:] + U[:q])
        prim = U[:p]
        G = least_rotation(prim)
        r = next(r for r in range(p) if prim[r:] + prim[:r] == G)
        L0 = len(start)
        tau = (-L0 - r) % p          # group coordinate = walk coordinate + tau
        # walk coordinates: T(i) = U[(i - L0) % P]; start window [0, L0)
        x, e, Z = 0, L0, [None] * P
        for k, _ in steps:
            if k == "u":
                e += 1
            else:
                Z[x] = e
                x += 1
        key = (p, G)
        if key not in groups_by_key:
            groups_by_key[key] = len(groups)
            groups.append(Group(p, G))
        gid = groups_by_key[key]
        gwalks[gid].append((P, Z, tau))

    for gid, g in enumerate(groups):
        p = g.p
        g.t = max(range(p), key=lambda i: g.G[i:] + g.G[:i])
        cols = []
        for x in range(p):
            Ex = []
            for P, Z, tau in gwalks[gid]:
                for rr in range(P // p):
                    xg = x + rr * p              # group coord
                    xw = xg - tau                # walk coord
                    zw = Z[xw % P] + P * (xw // P)
                    Ex.append(zw + tau - rr * p)
            cols.append(sorted(Ex))
        h = len(cols[0])
        for i in range(h):
            layers.append(Layer(gid, [cols[x][i] for x in range(p)]))
            g.layers.append(len(layers) - 1)
    return groups, layers


# --------------------------------------------------------------------------
# Sections 4-6: connecting everything to epsilon within a second budget W
# --------------------------------------------------------------------------
class NeedRule(Exception):
    """Section 5.3 found no record: the rule (2.6) it relies on is missing from m."""


@dataclass
class Record:
    target: int   # layer id D
    aD: int       # window [aD, bD) in D's group coordinates
    bD: int
    delta: int    # source coordinate = D coordinate + delta


class Connector:
    def __init__(self, groups, layers, walks, V, stats):
        self.G, self.L, self.W, self.V, self.stats = groups, layers, walks, V, stats
        self.rooted = set()
        self.used_rules = []
        self.requests = defaultdict(list)   # host layer -> [(child gid, record, child layer, e child coords)]
        self.blocks = []

    # layer helpers --------------------------------------------------------
    def p(self, lid):
        return self.G[self.L[lid].gid].p

    def f(self, lid, x):
        return self.L[lid].z(x - 1, self.p(lid))

    def l(self, lid, x):
        return self.L[lid].z(x, self.p(lid))

    def min_vertex(self, lid):
        g = self.G[self.L[lid].gid]
        x = min(range(g.p), key=lambda x: self.f(lid, x) - x)
        return x, self.f(lid, x)

    def find_record(self, gid, t, end):
        """Section 5.3 / 7: an actual window on some layer, agreeing with A throughout,
        containing w = A[t:end] at an alignment different from A's own."""
        g = self.G[gid]
        w = g.word(t, end)
        for D, lay in enumerate(self.L):
            gd = self.G[lay.gid]
            for x in range(gd.p):
                lo, hi = self.f(D, x), self.l(D, x)
                for e in range(max(lo, x + len(w)), hi + 1):
                    word = gd.word(x, e)
                    o = word.find(w)
                    while o != -1:
                        delta = t - (x + o)
                        same = lay.gid == gid and delta % g.p == 0
                        if not same and all(gd.A(y) == g.A(y + delta) for y in range(x, e)):
                            return Record(D, x, e, delta)
                        o = word.find(w, o + 1)
        raise AssertionError("no record found (Lemma 5.1 violated)")

    # Lemma 4.2 -------------------------------------------------------------
    def band_join(self, g, t, E, H, h):
        if E < H:
            assert H - E <= h
            self.W.follow(g.A, [(t, E), (t, H), (t + h, E + h)])
        self.stats["band_join"] += 1

    # Lemma 4.4 -------------------------------------------------------------
    def fulfill(self, host):
        gh = self.G[self.L[host].gid]
        q = gh.p
        by_child = defaultdict(list)
        for cg, rec, child, e in self.requests[host]:
            by_child[cg].append((rec, child, e))
        recs = []  # [cg, p_i, a, b, delta(child = host + delta), sorted ends(host coords), children]
        for cg in sorted(by_child, key=lambda c: (self.G[c].p, self.G[c].G)):
            items = by_child[cg]
            rec = items[0][0]
            assert all(it[0] is rec for it in items), "child group used two records"
            ends = sorted(e - rec.delta for _, _, e in items)
            recs.append([cg, self.G[cg].p, rec.aD, rec.bD, rec.delta, ends, [c for _, c, _ in items]])
        for cg, pi, a, b, dl, ends, _ in recs:
            for e in ends:
                assert a + pi < e <= b - pi, "(4.2) violated"
        j = min(range(len(recs)), key=lambda i: recs[i][1])
        rj = recs[j]
        aj, bj = rj[2], rj[3]
        others = []
        for i, r in enumerate(recs):
            if i == j:
                continue
            # translate by a multiple of q to lie between record j and its next copy
            s = ((aj - r[2]) // q) * q
            while (r[2] + s, r[3] + s) < (aj, bj):
                s += q
            while (r[2] + s, r[3] + s) >= (aj + q, bj + q):
                s -= q
            r[2] += s; r[3] += s; r[4] -= s
            r[5] = [e + s for e in r[5]]
            others.append(r)
        others.sort(key=lambda r: (r[2], r[3]))
        T = aj + q
        lst = [(aj, e) for e in rj[5]] + [(aj, bj)]
        targets = []
        for r in others:
            cg, pi, ai, bi, dl, ends, _ = r
            for e in ends:
                Q = (min(ai + pi, T), max(e, bj))
                lst.append(Q)
                targets.append((r, Q, e))
        last = lst[-1]
        lst.append((aj + q, rj[5][0] + q))
        self.W.follow(gh.A, lst)
        contact = gh.word(T, last[1]) if T <= last[1] else ""
        assert len(contact) < rj[1]
        self.W.round_trip(contact)
        for r, Q, e in targets:
            cg, pi, ai, bi, dl, _, _ = r
            gc = self.G[cg]
            Ai = (lambda gc, dl: (lambda y: gc.A(y + dl)))(gc, dl)
            self.W.follow(Ai, [Q, (ai + pi, e + pi), (Q[0] + pi, Q[1] + pi)])
        self.rooted.add(host)
        for r in recs:
            self.rooted.update(r[6])
        self.stats["host_requests"] += 1

    # Section 5 --------------------------------------------------------------
    def process(self):
        order = sorted(range(len(self.G)), key=lambda i: (self.G[i].p, self.G[i].G))
        for gid in order:
            self.process_group(gid)
        self.open_cycles()

    def process_group(self, gid):
        g = self.G[gid]
        p, t, Ls = g.p, g.t, g.layers
        recipients = [i for i, lid in enumerate(Ls) if self.requests[lid]]
        for i in recipients:
            self.fulfill(Ls[i])
        upper = Ls
        if recipients:
            bi = max(recipients)
            base = Ls[bi]
            _, rec, _, _ = self.requests[base][0]
            a, h = rec.aD, rec.bD           # baseline vertex S of length < 2p
            assert h - a < 2 * p
            for lid in Ls[:bi]:
                if lid in self.rooted:
                    continue
                x, fx = self.min_vertex(lid)
                if fx - x <= p:
                    self.W.round_trip(g.word(x, fx))
                else:
                    e = min(self.l(lid, a), h)
                    assert p < e - a
                    self.W.follow(g.A, [(a, h), (a + p, e + p), (a + p, h + p)])
                self.rooted.add(lid)
                self.stats["below_baseline"] += 1
            upper = Ls[bi + 1:]
        n = len(upper)
        if n == 0:
            return
        top = upper[-1]
        H = self.f(top, t)
        # 5.2 easy connections
        if H - t <= n * p:
            self.W.round_trip(g.word(t, H))
            self.rooted.update(upper)
            self.stats["easy_round_trip"] += 1
            return
        if recipients and self.l(base, t) >= H - n * p:
            self.band_join(g, t, self.l(base, t), H, n * p)
            self.rooted.update(upper)
            return
        # 5.3 threshold search
        k = 1
        while sum(1 for lid in upper if self.l(lid, t) >= H - k * p) > k:
            k += 1
        K = k * p
        topk = upper[-k:]
        Lstar = topk[0]
        assert self.l(Lstar, t) >= H - (k - 1) * p
        self.used_rules.append((g.word(t, H), p, k))      # rule (2.6) relied on here
        try:
            rec = self.find_record(gid, t, H - K)
        except AssertionError:
            raise NeedRule((g.word(t, H), p, k))
        a, b = rec.aD + rec.delta, rec.bD + rec.delta     # record in A coordinates
        q = self.p(rec.target)
        assert a <= t and b >= H - K and b - a < p + q
        # 5.4 distributing budgets
        rH, rL = self.f(top, a), self.l(Lstar, a)
        assert rH - rL <= K and rH <= b + K
        rest = upper
        if rL + K >= b:
            self.blocks.append(dict(kind="collective", gid=gid, layers=topk, rec=rec, k=k, t=t, H=H))
            rest = upper[:-k]
        for lid in rest:
            x, fx = self.min_vertex(lid)
            if fx - x <= p:
                self.W.round_trip(g.word(x, fx))
                self.rooted.add(lid)
                self.stats["upper_round_trip"] += 1
            elif self.l(lid, a) + p >= b:
                self.blocks.append(dict(kind="individual", gid=gid, layers=[lid], rec=rec, k=1, t=t))
            else:
                e = self.f(lid, a)
                assert a + p < e <= self.l(lid, a) < b - p
                assert q > p
                self.requests[rec.target].append((gid, rec, lid, e))
                self.stats["requests"] += 1

    # link realisations -------------------------------------------------------
    def link(self, blk, short_q=None):
        """Execute a block's planned link. With short_q, realise it so the piece has a
        vertex of length <= short_q (Lemma 6.2) and return that contact word."""
        g = self.G[blk["gid"]]
        p, rec = g.p, blk["rec"]
        a, b = rec.aD + rec.delta, rec.bD + rec.delta
        if blk["kind"] == "individual":
            lid = blk["layers"][0]
            e2 = max(b, self.f(lid, a + p))
            assert e2 <= min(b + p, self.l(lid, a + p))
            self.W.follow(g.A, [(a, b), (a + p, e2), (a + p, b + p)])
            self.stats["individual_link"] += 1
            return g.word(a + p, b) if a + p <= b else ""
        K = blk["k"] * p
        lids = blk["layers"]                       # ascending; top = last
        if short_q is not None:
            for lid in lids:
                x, fx = self.min_vertex(lid)
                if fx - x <= short_q:
                    self._collective(g, lids, a, b, K)
                    return g.word(x, fx)
            t = blk["t"]
            x = min(t, a + K)
            F = self.f(lids[-1], x)
            assert b - x < short_q and b < F <= blk["H"] <= b + K
            self.W.follow(g.A, [(a, b), (x, b), (x, F), (a + K, b + K)])
            self.stats["short_contact_link"] += 1
            return g.word(x, b)
        self._collective(g, lids, a, b, K)
        return None

    def _collective(self, g, lids, a, b, K):
        """Lemma 4.3."""
        Y = max(b, self.f(lids[-1], a))
        wins = [(a, b)]
        for lid in reversed(lids):
            x = next(x for x in range(a, a + K + 1) if self.l(lid, x) >= Y)
            wins.append((x, Y))
        wins.append((a + K, b + K))
        self.W.follow(g.A, wins)
        self.stats["collective_link"] += 1

    def internal_join(self, blk):
        """Lemma 6.1."""
        if blk["kind"] == "collective" and blk["k"] > 1:
            g = self.G[blk["gid"]]
            t, lids = blk["t"], blk["layers"]
            self.band_join(g, t, self.l(lids[0], t), self.f(lids[-1], t), (blk["k"] - 1) * g.p)

    # Section 6 ---------------------------------------------------------------
    def open_cycles(self):
        block_of = {}
        for bi, blk in enumerate(self.blocks):
            for lid in blk["layers"]:
                assert lid not in self.rooted and lid not in block_of
                block_of[lid] = bi
        nxt = [block_of.get(blk["rec"].target) for blk in self.blocks]
        state = [0] * len(self.blocks)   # 0 new, 1 on stack, 2 done
        on_cycle, omitted, handled = set(), set(), set()
        maxp = max(g.p for g in self.G)
        rotkey = lambda bi: (lambda g: ((g.G[g.t:] + g.G[:g.t]) * (2 * maxp // g.p + 2))[: 2 * maxp])(self.G[self.blocks[bi]["gid"]])
        for s in range(len(self.blocks)):
            path, v = [], s
            while v is not None and state[v] == 0:
                state[v] = 1
                path.append(v)
                v = nxt[v]
            if v is not None and state[v] == 1:
                cyc = path[path.index(v):]
                on_cycle.update(cyc)
                if len(cyc) == 1:
                    blk = self.blocks[cyc[0]]
                    self.internal_join(blk)
                    rec = blk["rec"]
                    g = self.G[blk["gid"]]
                    R = g.word(rec.aD + rec.delta, rec.bD + rec.delta)
                    assert len(R) < g.p
                    self.W.round_trip(R)
                    omitted.add(cyc[0])
                    self.stats["self_cycle"] += 1
                else:
                    sel = max(cyc, key=rotkey)
                    tb = nxt[sel]
                    q = self.G[self.blocks[tb]["gid"]].p
                    contact = self.link(self.blocks[sel], short_q=q)
                    assert len(contact) <= q
                    self.internal_join(self.blocks[tb])
                    self.W.round_trip(contact)
                    omitted.add(tb)
                    handled.add(sel)
                    self.stats["opened_cycle"] += 1
            for u in path:
                state[u] = 2
        for bi, blk in enumerate(self.blocks):
            if bi not in omitted and bi not in handled:
                self.link(blk)


# --------------------------------------------------------------------------
# Euler tour and driver
# --------------------------------------------------------------------------
def euler_superstring(up, down):
    out = defaultdict(list)
    for s in sorted(up):
        out[s[:-1]].append([s, s[-1], up[s]])
    for s in sorted(down):
        out[s].append([s[1:], None, down[s]])
    ptr = defaultdict(int)
    stack, letters = [("", None)], []
    while stack:
        v, _ = stack[-1]
        lst = out[v]
        while ptr[v] < len(lst) and lst[ptr[v]][2] == 0:
            ptr[v] += 1
        if ptr[v] < len(lst):
            e = lst[ptr[v]]
            e[2] -= 1
            stack.append((e[0], e[1]))
        else:
            letters.append(stack.pop()[1])
    return "".join(ch for ch in reversed(letters) if ch)


def connected_to_eps(up, down):
    adj = defaultdict(set)
    for s in up:
        adj[s[:-1]].add(s); adj[s].add(s[:-1])
    for s in down:
        adj[s].add(s[1:]); adj[s[1:]].add(s)
    seen, st = {""}, [""]
    while st:
        for w in adj[st.pop()]:
            if w not in seen:
                seen.add(w); st.append(w)
    return all(v in seen for v in adj)


@dataclass
class Result:
    superstring: str
    W: int
    stats: dict


def two_approx(strings, rule_filter=None):
    req = reduce_inputs(strings)
    if not req:
        return Result("", 0, {})
    C = Counts(req, rule_filter)
    up, down = C.edges()
    walks = closed_walks(up, down)
    groups, layers = build_layers(walks)
    if CHECK:
        for lay in layers:
            g = groups[lay.gid]
            for x in range(g.p):
                for e in range(lay.z(x - 1, g.p), lay.z(x, g.p) + 1):
                    assert g.word(x, e) in C.V
    stats = Counter()
    added = Walks(C.V)
    conn = Connector(groups, layers, added, C.V, stats)
    conn.process()
    assert added.cost <= C.W, ("budget exceeded", added.cost, C.W)
    U, D = up + added.up, down + added.down
    assert connected_to_eps(U, D), "result not connected"
    T = euler_superstring(U, D)
    assert len(T) <= 2 * C.W
    if CHECK:
        assert all(s in T for s in strings if s)
    stats.update(groups=len(groups), layers=len(layers), added_cost=added.cost)
    res = Result(T, C.W, dict(stats))
    res.used_rules, res.counts = conn.used_rules, C
    return res


if __name__ == "__main__":
    for inst in [["cde", "abc", "efg"], ["ab", "ba", "bb"], ["aaaa"], ["abab", "baba", "bbb"]]:
        r = two_approx(inst)
        print(inst, "->", r.superstring, "W =", r.W, r.stats)


def two_approx_lazy(strings, max_rounds=10_000, batch=True):
    """Same guarantee (len <= 2W' <= 2 OPT) using only extension rules plus the
    periodic rules that Section 5.3 actually needs, added on demand.

    Any subset of rules (2.6) still gives a lower bound W' <= OPT (Lemma 2.1), and
    the connection argument uses a periodic rule only in Section 5.3."""
    rules, texts = set(), set()   # texts: (period, least rotation) enabled wholesale
    keep = lambda v, p, k: (v, p, k) in rules or (p, least_rotation(v[:p])) in texts
    for rnd in range(max_rounds):
        try:
            r = two_approx(strings, rule_filter=keep)
            r.stats["lazy_rules"] = len(rules) + len(texts)
            r.stats["lazy_rounds"] = rnd + 1
            return r
        except NeedRule as e:
            v, p, k = e.args[0]
            assert not keep(v, p, k), "rule already present but still no record"
            if batch:
                texts.add((p, least_rotation(v[:p])))
            else:
                rules.add((v, p, k))
    raise RuntimeError("lazy rule loop did not converge")
