//! Sections 3-7 of openai/math #128 on (state, length) vertices (port of scs2.py).

use crate::counts::Counts;
use crate::sam::{key, kl, kx, FxMap, Sam, EPS};
use std::collections::BTreeMap;

#[inline]
fn md(a: i64, b: i64) -> i64 {
    a.rem_euclid(b)
}
#[inline]
fn dv(a: i64, b: i64) -> i64 {
    a.div_euclid(b)
}

pub struct NeedRule(pub usize, pub Vec<u8>); // (p, least rotation) to enable

// ---------------------------------------------------------------- string utils
fn least_period(u: &[u8]) -> usize {
    let n = u.len();
    let mut pi = vec![0usize; n];
    for i in 1..n {
        let mut k = pi[i - 1];
        while k > 0 && u[i] != u[k] {
            k = pi[k - 1];
        }
        if u[i] == u[k] {
            k += 1;
        }
        pi[i] = k;
    }
    let c = n - pi[n - 1];
    if n % c == 0 { c } else { n }
}

/// Start index of the least (or greatest) rotation (two-pointer minimum expression).
fn extreme_rotation(s: &[u8], greatest: bool) -> usize {
    let n = s.len();
    let (mut i, mut j, mut k) = (0usize, 1usize, 0usize);
    while i < n && j < n && k < n {
        let (x, y) = (s[(i + k) % n], s[(j + k) % n]);
        if x == y {
            k += 1;
            continue;
        }
        if (x > y) != greatest { i += k + 1 } else { j += k + 1 }
        if i == j {
            j += 1;
        }
        k = 0;
    }
    i.min(j)
}

// ---------------------------------------------------------------- groups / layers
pub struct Group {
    pub p: i64,
    pub g: Vec<u8>,
    pub t: i64,
    pub layers: Vec<usize>,
}
impl Group {
    #[inline]
    fn a(&self, i: i64) -> u8 {
        self.g[md(i, self.p) as usize]
    }
}

pub struct Layer {
    pub gid: usize,
    pub zb: Vec<i64>,
}

#[derive(Clone, Copy)]
struct Text {
    gid: usize,
    delta: i64,
}

#[derive(Clone, Copy, PartialEq)]
struct Record {
    target: usize,
    ad: i64,
    bd: i64,
    delta: i64,
}

#[derive(PartialEq, Clone, Copy)]
enum Kind {
    Collective,
    Individual,
}
struct Block {
    kind: Kind,
    gid: usize,
    layers: Vec<usize>,
    rec: Record,
    k: i64,
    t: i64,
    h: i64,
}
struct Request {
    cg: usize,
    rec: Record,
    child: usize,
    e: i64,
}

pub struct Ctx<'a> {
    pub sam: &'a Sam,
    pub groups: Vec<Group>,
    pub layers: Vec<Layer>,
    pub up: FxMap<u64, u32>,
    pub down: FxMap<u64, u32>,
    pub cost: i64,
    rooted: Vec<bool>,
    requests: Vec<Vec<Request>>,
    blocks: Vec<Block>,
    pub stats: BTreeMap<&'static str, u64>,
}

impl<'a> Ctx<'a> {
    #[inline]
    fn ch(&self, t: Text, i: i64) -> u8 {
        self.groups[t.gid].a(i + t.delta)
    }
    fn p(&self, lid: usize) -> i64 {
        self.groups[self.layers[lid].gid].p
    }
    fn z(&self, lid: usize, x: i64) -> i64 {
        let p = self.p(lid);
        self.layers[lid].zb[md(x, p) as usize] + p * dv(x, p)
    }
    fn f(&self, lid: usize, x: i64) -> i64 {
        self.z(lid, x - 1)
    }
    fn l(&self, lid: usize, x: i64) -> i64 {
        self.z(lid, x)
    }
    fn min_vertex(&self, lid: usize) -> (i64, i64) {
        let p = self.p(lid);
        let x = (0..p).min_by_key(|&x| self.f(lid, x) - x).unwrap();
        (x, self.f(lid, x))
    }
    fn win(&self, t: Text, a: i64, b: i64) -> u64 {
        self.sam
            .vertex((a..b).map(|i| self.ch(t, i)))
            .unwrap_or_else(|| panic!("window not a vertex"))
    }
    fn bump(&mut self, up: bool, v: u64) {
        let m = if up { &mut self.up } else { &mut self.down };
        *m.entry(v).or_insert(0) += 1;
    }

    /// Lemma 4.1: follow ordered windows of one text, recording edges.
    fn follow(&mut self, t: Text, wins: &[(i64, i64)]) {
        let first = self.win(t, wins[0].0, wins[0].1);
        let last = self.win(t, wins[wins.len() - 1].0, wins[wins.len() - 1].1);
        assert_eq!(first, last, "walk not closed");
        for w in wins.windows(2) {
            let ((a, b), (c, d)) = (w[0], w[1]);
            assert!(a <= c && b <= d, "unordered windows");
            let mut v = self.win(t, a, b);
            if c <= b {
                for _ in a..c {
                    self.bump(false, v);
                    v = self.sam.suf(v);
                }
            } else {
                for _ in a..b {
                    self.bump(false, v);
                    v = self.sam.suf(v);
                }
                for i in b..c {
                    let u = self.sam.ext(EPS, self.ch(t, i)).unwrap();
                    self.bump(true, u);
                    self.bump(false, u);
                }
                v = EPS;
            }
            for e in b.max(c) + 1..=d {
                v = self.sam.ext(v, self.ch(t, e - 1)).expect("window not a vertex");
                self.bump(true, v);
            }
            self.cost += d - b;
        }
    }

    /// epsilon -> word(a, b) -> epsilon, cost b - a.
    fn round_trip(&mut self, t: Text, a: i64, b: i64) {
        let mut v = EPS;
        for i in a..b {
            v = self.sam.ext(v, self.ch(t, i)).expect("round trip word not a vertex");
            self.bump(true, v);
        }
        for _ in a..b {
            self.bump(false, v);
            v = self.sam.suf(v);
        }
        self.cost += b - a;
    }

    fn gt(gid: usize) -> Text {
        Text { gid, delta: 0 }
    }
    fn st(&mut self, k: &'static str) {
        *self.stats.entry(k).or_insert(0) += 1;
    }

    // ------------------------------------------------------------ Section 5.3 search
    fn find_record(&self, gid: usize, t: i64, end: i64) -> Option<Record> {
        let g = &self.groups[gid];
        let w: Vec<u8> = (t..end).map(|i| g.a(i)).collect();
        let wl = w.len() as i64;
        // Any record satisfying the conditions is valid (Section 5.3 / 7); the scan order
        // can be varied via SCS_RECORD_ORDER=rev|same to exercise different branches.
        let mut ids: Vec<usize> = (0..self.layers.len()).collect();
        match std::env::var("SCS_RECORD_ORDER").as_deref() {
            Ok("rev") => ids.reverse(),
            Ok("same") => ids.sort_by_key(|&d| (self.layers[d].gid != gid, std::cmp::Reverse(d))),
            _ => {}
        }
        for d in ids {
            let gd = &self.groups[self.layers[d].gid];
            for x in 0..gd.p {
                let (lo, hi) = (self.f(d, x), self.l(d, x));
                if hi - x < wl {
                    continue;
                }
                let word: Vec<u8> = (x..hi).map(|i| gd.a(i)).collect();
                for e in lo.max(x + wl)..=hi {
                    let wn = (e - x) as usize;
                    for o in 0..=(wn - w.len()) {
                        if word[o..o + w.len()] != w[..] {
                            continue;
                        }
                        let delta = t - (x + o as i64);
                        let same = self.layers[d].gid == gid && md(delta, g.p) == 0;
                        if !same && (x..e).all(|y| gd.a(y) == g.a(y + delta)) {
                            return Some(Record { target: d, ad: x, bd: e, delta });
                        }
                    }
                }
            }
        }
        None
    }

    fn band_join(&mut self, gid: usize, t: i64, e: i64, h: i64, budget: i64) {
        if e < h {
            assert!(h - e <= budget);
            self.follow(Self::gt(gid), &[(t, e), (t, h), (t + budget, e + budget)]);
        }
        self.st("band_join");
    }

    // ------------------------------------------------------------ Lemma 4.4
    fn fulfill(&mut self, host: usize) {
        let hg = self.layers[host].gid;
        let q = self.groups[hg].p;
        let reqs = std::mem::take(&mut self.requests[host]);
        // child groups ordered by (p, text)
        let mut cgs: Vec<usize> = reqs.iter().map(|r| r.cg).collect();
        cgs.sort_by(|&a, &b| (self.groups[a].p, &self.groups[a].g).cmp(&(self.groups[b].p, &self.groups[b].g)));
        cgs.dedup();
        struct R {
            pi: i64,
            a: i64,
            b: i64,
            dl: i64,
            cg: usize,
            ends: Vec<i64>,
            children: Vec<usize>,
        }
        let mut recs: Vec<R> = Vec::new();
        for &cg in &cgs {
            let items: Vec<&Request> = reqs.iter().filter(|r| r.cg == cg).collect();
            let rec = items[0].rec;
            assert!(items.iter().all(|r| r.rec == rec), "child group used two records");
            let mut ends: Vec<i64> = items.iter().map(|r| r.e - rec.delta).collect();
            ends.sort_unstable();
            recs.push(R { pi: self.groups[cg].p, a: rec.ad, b: rec.bd, dl: rec.delta, cg, ends,
                          children: items.iter().map(|r| r.child).collect() });
        }
        for r in &recs {
            for &e in &r.ends {
                assert!(r.a + r.pi < e && e <= r.b - r.pi, "(4.2) violated");
            }
        }
        let j = (0..recs.len()).min_by_key(|&i| recs[i].pi).unwrap();
        let (aj, bj, pj) = (recs[j].a, recs[j].b, recs[j].pi);
        let mut others: Vec<usize> = Vec::new();
        for i in 0..recs.len() {
            if i == j {
                continue;
            }
            let r = &mut recs[i];
            let mut s = dv(aj - r.a, q) * q;
            while (r.a + s, r.b + s) < (aj, bj) {
                s += q;
            }
            while (r.a + s, r.b + s) >= (aj + q, bj + q) {
                s -= q;
            }
            r.a += s;
            r.b += s;
            r.dl -= s;
            for e in r.ends.iter_mut() {
                *e += s;
            }
            others.push(i);
        }
        others.sort_by_key(|&i| (recs[i].a, recs[i].b));
        let tt = aj + q;
        let mut lst: Vec<(i64, i64)> = recs[j].ends.iter().map(|&e| (aj, e)).collect();
        lst.push((aj, bj));
        let mut targets: Vec<(usize, (i64, i64), i64)> = Vec::new();
        for &i in &others {
            let r = &recs[i];
            for &e in &r.ends {
                let qq = ((r.a + r.pi).min(tt), e.max(bj));
                lst.push(qq);
                targets.push((i, qq, e));
            }
        }
        let last = *lst.last().unwrap();
        lst.push((aj + q, recs[j].ends[0] + q));
        self.follow(Self::gt(hg), &lst);
        if tt <= last.1 {
            assert!(last.1 - tt < pj);
            self.round_trip(Self::gt(hg), tt, last.1);
        }
        for (i, qq, e) in targets {
            let r = &recs[i];
            let ti = Text { gid: r.cg, delta: r.dl };
            let pi = r.pi;
            let a = r.a;
            self.follow(ti, &[qq, (a + pi, e + pi), (qq.0 + pi, qq.1 + pi)]);
        }
        self.rooted[host] = true;
        for r in &recs {
            for &c in &r.children {
                self.rooted[c] = true;
            }
        }
        self.st("host_requests");
    }

    // ------------------------------------------------------------ Section 5
    fn process_group(&mut self, gid: usize) -> Result<(), NeedRule> {
        let (p, t) = (self.groups[gid].p, self.groups[gid].t);
        let ls = self.groups[gid].layers.clone();
        let g = Self::gt(gid);
        let recipients: Vec<usize> = (0..ls.len()).filter(|&i| !self.requests[ls[i]].is_empty()).collect();
        // baseline record must be read before fulfill() consumes the requests
        let base_rec = recipients.last().map(|&bi| self.requests[ls[bi]][0].rec);
        for &i in &recipients {
            self.fulfill(ls[i]);
        }
        let mut upper: Vec<usize> = ls.clone();
        let mut base: Option<usize> = None;
        if let Some(&bi) = recipients.last() {
            base = Some(ls[bi]);
            let rec = base_rec.unwrap();
            let (a, h) = (rec.ad, rec.bd);
            assert!(h - a < 2 * p);
            for &lid in &ls[..bi] {
                if self.rooted[lid] {
                    continue;
                }
                let (x, fx) = self.min_vertex(lid);
                if fx - x <= p {
                    self.round_trip(g, x, fx);
                } else {
                    let e = self.l(lid, a).min(h);
                    assert!(p < e - a);
                    self.follow(g, &[(a, h), (a + p, e + p), (a + p, h + p)]);
                }
                self.rooted[lid] = true;
                self.st("below_baseline");
            }
            upper = ls[bi + 1..].to_vec();
        }
        let n = upper.len() as i64;
        if n == 0 {
            return Ok(());
        }
        let top = *upper.last().unwrap();
        let hh = self.f(top, t);
        if hh - t <= n * p {
            self.round_trip(g, t, hh);
            for &lid in &upper {
                self.rooted[lid] = true;
            }
            self.st("easy_round_trip");
            return Ok(());
        }
        if let Some(b) = base {
            if self.l(b, t) >= hh - n * p {
                let e = self.l(b, t);
                self.band_join(gid, t, e, hh, n * p);
                for &lid in &upper {
                    self.rooted[lid] = true;
                }
                return Ok(());
            }
        }
        let mut k = 1i64;
        while upper.iter().filter(|&&lid| self.l(lid, t) >= hh - k * p).count() as i64 > k {
            k += 1;
        }
        let kk = k * p;
        let topk: Vec<usize> = upper[(n - k) as usize..].to_vec();
        let lstar = topk[0];
        assert!(self.l(lstar, t) >= hh - (k - 1) * p);
        let rec = match self.find_record(gid, t, hh - kk) {
            Some(r) => r,
            None => return Err(NeedRule(p as usize, self.groups[gid].g.clone())),
        };
        let (a, b) = (rec.ad + rec.delta, rec.bd + rec.delta);
        let q = self.p(rec.target);
        assert!(a <= t && b >= hh - kk && b - a < p + q);
        let (rh, rl) = (self.f(top, a), self.l(lstar, a));
        assert!(rh - rl <= kk && rh <= b + kk);
        let mut rest = upper.clone();
        if rl + kk >= b {
            self.blocks.push(Block { kind: Kind::Collective, gid, layers: topk, rec, k, t, h: hh });
            rest = upper[..(n - k) as usize].to_vec();
        }
        for lid in rest {
            let (x, fx) = self.min_vertex(lid);
            if fx - x <= p {
                self.round_trip(g, x, fx);
                self.rooted[lid] = true;
                self.st("upper_round_trip");
            } else if self.l(lid, a) + p >= b {
                self.blocks.push(Block { kind: Kind::Individual, gid, layers: vec![lid], rec, k: 1, t, h: hh });
            } else {
                let e = self.f(lid, a);
                assert!(a + p < e && e <= self.l(lid, a) && self.l(lid, a) < b - p);
                assert!(q > p);
                self.requests[rec.target].push(Request { cg: gid, rec, child: lid, e });
                self.st("requests");
            }
        }
        Ok(())
    }

    fn collective(&mut self, gid: usize, lids: &[usize], a: i64, b: i64, kk: i64) {
        let y = b.max(self.f(*lids.last().unwrap(), a));
        let mut wins = vec![(a, b)];
        for &lid in lids.iter().rev() {
            let x = (a..=a + kk).find(|&x| self.l(lid, x) >= y).unwrap();
            wins.push((x, y));
        }
        wins.push((a + kk, b + kk));
        self.follow(Self::gt(gid), &wins);
        self.st("collective_link");
    }

    /// Execute a planned link; with short_q, realise Lemma 6.2 and return the contact window.
    fn link(&mut self, bi: usize, short_q: Option<i64>) -> Option<(i64, i64)> {
        let (kind, gid, rec, k, t, hh) = {
            let b = &self.blocks[bi];
            (b.kind, b.gid, b.rec, b.k, b.t, b.h)
        };
        let lids = self.blocks[bi].layers.clone();
        let p = self.groups[gid].p;
        let g = Self::gt(gid);
        let (a, b) = (rec.ad + rec.delta, rec.bd + rec.delta);
        if kind == Kind::Individual {
            let lid = lids[0];
            let e2 = b.max(self.f(lid, a + p));
            assert!(e2 <= (b + p).min(self.l(lid, a + p)));
            self.follow(g, &[(a, b), (a + p, e2), (a + p, b + p)]);
            self.st("individual_link");
            return Some((a + p, (a + p).max(b)));
        }
        let kk = k * p;
        if let Some(sq) = short_q {
            for &lid in &lids {
                let (x, fx) = self.min_vertex(lid);
                if fx - x <= sq {
                    self.collective(gid, &lids, a, b, kk);
                    return Some((x, fx));
                }
            }
            let x = t.min(a + kk);
            let ff = self.f(*lids.last().unwrap(), x);
            assert!(b - x < sq && b < ff && ff <= hh && hh <= b + kk);
            self.follow(g, &[(a, b), (x, b), (x, ff), (a + kk, b + kk)]);
            self.st("short_contact_link");
            return Some((x, b));
        }
        self.collective(gid, &lids, a, b, kk);
        None
    }

    fn internal_join(&mut self, bi: usize) {
        let b = &self.blocks[bi];
        if b.kind == Kind::Collective && b.k > 1 {
            let (gid, t, k) = (b.gid, b.t, b.k);
            let (lo, hi) = (b.layers[0], *b.layers.last().unwrap());
            let p = self.groups[gid].p;
            let (e, h) = (self.l(lo, t), self.f(hi, t));
            self.band_join(gid, t, e, h, (k - 1) * p);
        }
    }

    /// Compare distinguished forward rotations of two groups (infinite words).
    fn rot_cmp(&self, g1: usize, g2: usize) -> std::cmp::Ordering {
        let (a, b) = (&self.groups[g1], &self.groups[g2]);
        let n = a.p + b.p;
        for i in 0..n {
            let (x, y) = (a.a(a.t + i), b.a(b.t + i));
            if x != y {
                return x.cmp(&y);
            }
        }
        std::cmp::Ordering::Equal
    }

    // ------------------------------------------------------------ Section 6
    fn open_cycles(&mut self) {
        let nb = self.blocks.len();
        let mut block_of: FxMap<usize, usize> = FxMap::default();
        for (bi, b) in self.blocks.iter().enumerate() {
            for &lid in &b.layers {
                assert!(!self.rooted[lid] && !block_of.contains_key(&lid));
                block_of.insert(lid, bi);
            }
        }
        let nxt: Vec<Option<usize>> = self.blocks.iter().map(|b| block_of.get(&b.rec.target).copied()).collect();
        let mut state = vec![0u8; nb];
        let mut omitted = vec![false; nb];
        let mut handled = vec![false; nb];
        for s in 0..nb {
            let mut path = Vec::new();
            let mut v = Some(s);
            while let Some(u) = v {
                if state[u] != 0 {
                    break;
                }
                state[u] = 1;
                path.push(u);
                v = nxt[u];
            }
            if let Some(u) = v {
                if state[u] == 1 {
                    let cyc: Vec<usize> = path[path.iter().position(|&x| x == u).unwrap()..].to_vec();
                    if cyc.len() == 1 {
                        let c = cyc[0];
                        self.internal_join(c);
                        let (gid, rec) = (self.blocks[c].gid, self.blocks[c].rec);
                        let (a, b) = (rec.ad + rec.delta, rec.bd + rec.delta);
                        assert!(b - a < self.groups[gid].p);
                        self.round_trip(Self::gt(gid), a, b);
                        omitted[c] = true;
                        self.st("self_cycle");
                    } else {
                        let mut sel = cyc[0];
                        for &c in &cyc[1..] {
                            if self.rot_cmp(self.blocks[c].gid, self.blocks[sel].gid) == std::cmp::Ordering::Greater {
                                sel = c;
                            }
                        }
                        let tb = nxt[sel].unwrap();
                        let q = self.groups[self.blocks[tb].gid].p;
                        let (ca, cb) = self.link(sel, Some(q)).unwrap();
                        assert!(cb - ca <= q);
                        self.internal_join(tb);
                        let sg = self.blocks[sel].gid;
                        self.round_trip(Self::gt(sg), ca, cb);
                        omitted[tb] = true;
                        handled[sel] = true;
                        self.st("opened_cycle");
                    }
                }
            }
            for u in path {
                state[u] = 2;
            }
        }
        for bi in 0..nb {
            if !omitted[bi] && !handled[bi] {
                self.link(bi, None);
            }
        }
    }
}

// ---------------------------------------------------------------- driver pieces
type Steps = Vec<Option<u8>>; // Some(c): up appending c; None: down

fn take(sam: &Sam, up: &mut FxMap<u64, u32>, down: &mut FxMap<u64, u32>, v: u64) -> Option<(u64, Option<u8>)> {
    let (x, l) = (kx(v), kl(v));
    for &(c, y) in &sam.next[x as usize] {
        let k = key(y, l + 1);
        if let Some(m) = up.get_mut(&k) {
            if *m > 0 {
                *m -= 1;
                return Some((k, Some(c)));
            }
        }
    }
    if let Some(m) = down.get_mut(&v) {
        if *m > 0 {
            *m -= 1;
            return Some((sam.suf(v), None));
        }
    }
    None
}

fn closed_walks(sam: &Sam, up: &FxMap<u64, u32>, down: &FxMap<u64, u32>) -> Vec<(u32, Steps)> {
    let (mut ru, mut rd) = (up.clone(), down.clone());
    let mut starts: Vec<u64> = up.keys().chain(down.keys()).copied().collect();
    starts.push(EPS);
    starts.sort_unstable_by_key(|&k| (kl(k), k));
    starts.dedup();
    let mut walks = Vec::new();
    for &s in &starts {
        while let Some((mut v, c)) = take(sam, &mut ru, &mut rd, s) {
            let mut steps = vec![c];
            while v != s {
                let (nv, c) = take(sam, &mut ru, &mut rd, v).expect("graph not balanced");
                steps.push(c);
                v = nv;
            }
            walks.push((kl(s), steps));
        }
    }
    walks
}

fn build_layers(walks: &[(u32, Steps)]) -> (Vec<Group>, Vec<Layer>) {
    let mut gmap: FxMap<(usize, Vec<u8>), usize> = FxMap::default();
    let mut groups: Vec<Group> = Vec::new();
    let mut gw: Vec<Vec<(i64, Vec<i64>, i64)>> = Vec::new();
    for (l0, steps) in walks {
        let u: Vec<u8> = steps.iter().filter_map(|c| *c).collect();
        let pp = u.len() as i64;
        let p = least_period(&u);
        let prim = &u[..p];
        let r = extreme_rotation(prim, false);
        let g: Vec<u8> = prim[r..].iter().chain(prim[..r].iter()).copied().collect();
        let tau = md(-(*l0 as i64) - r as i64, p as i64);
        let mut z = vec![0i64; pp as usize];
        let (mut x, mut e) = (0usize, *l0 as i64);
        for c in steps {
            if c.is_some() {
                e += 1;
            } else {
                z[x] = e;
                x += 1;
            }
        }
        let gid = *gmap.entry((p, g.clone())).or_insert_with(|| {
            groups.push(Group { p: p as i64, g: g.clone(), t: 0, layers: Vec::new() });
            gw.push(Vec::new());
            groups.len() - 1
        });
        gw[gid].push((pp, z, tau));
    }
    let mut layers = Vec::new();
    for (gid, grp) in groups.iter_mut().enumerate() {
        let p = grp.p;
        grp.t = extreme_rotation(&grp.g, true) as i64;
        let mut cols: Vec<Vec<i64>> = Vec::with_capacity(p as usize);
        for x in 0..p {
            let mut ex = Vec::new();
            for (pp, z, tau) in &gw[gid] {
                for rr in 0..pp / p {
                    let xw = x + rr * p - tau;
                    let zw = z[md(xw, *pp) as usize] + pp * dv(xw, *pp);
                    ex.push(zw + tau - rr * p);
                }
            }
            ex.sort_unstable();
            cols.push(ex);
        }
        for i in 0..cols[0].len() {
            layers.push(Layer { gid, zb: cols.iter().map(|c| c[i]).collect() });
            grp.layers.push(layers.len() - 1);
        }
    }
    (groups, layers)
}

fn euler(sam: &Sam, up: &FxMap<u64, u32>, down: &FxMap<u64, u32>) -> Vec<u8> {
    let (mut ru, mut rd) = (up.clone(), down.clone());
    let mut stack: Vec<(u64, Option<u8>)> = vec![(EPS, None)];
    let mut out: Vec<Option<u8>> = Vec::new();
    while let Some(&(v, _)) = stack.last() {
        match take(sam, &mut ru, &mut rd, v) {
            Some((nv, c)) => stack.push((nv, c)),
            None => out.push(stack.pop().unwrap().1),
        }
    }
    assert!(ru.values().all(|&m| m == 0) && rd.values().all(|&m| m == 0), "result not connected");
    out.iter().rev().filter_map(|c| *c).collect()
}

pub struct Output {
    pub superstring: Vec<u8>,
    pub w: u64,
    pub stats: BTreeMap<&'static str, u64>,
}

/// One full run for fixed counts. Err(NeedRule) if Section 5.3 needs a periodic rule.
pub fn run(sam: &Sam, counts: &Counts, check: bool) -> Result<Output, NeedRule> {
    let (up, down) = counts.edges(sam);
    let walks = closed_walks(sam, &up, &down);
    let (groups, layers) = build_layers(&walks);
    let nl = layers.len();
    let mut cx = Ctx {
        sam,
        groups,
        layers,
        up: FxMap::default(),
        down: FxMap::default(),
        cost: 0,
        rooted: vec![false; nl],
        requests: (0..nl).map(|_| Vec::new()).collect(),
        blocks: Vec::new(),
        stats: BTreeMap::new(),
    };
    if check {
        for lid in 0..nl {
            let gid = cx.layers[lid].gid;
            for x in 0..cx.groups[gid].p {
                for e in cx.f(lid, x)..=cx.l(lid, x) {
                    cx.win(Ctx::gt(gid), x, e);
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..cx.groups.len()).collect();
    order.sort_by(|&a, &b| (cx.groups[a].p, &cx.groups[a].g).cmp(&(cx.groups[b].p, &cx.groups[b].g)));
    for gid in order {
        cx.process_group(gid)?;
    }
    cx.open_cycles();
    assert!(cx.cost as u64 <= counts.w, "budget exceeded");
    let (mut tu, mut td) = (up, down);
    for (k, v) in cx.up.drain() {
        *tu.entry(k).or_insert(0) += v;
    }
    for (k, v) in cx.down.drain() {
        *td.entry(k).or_insert(0) += v;
    }
    let t = euler(sam, &tu, &td);
    assert!(t.len() as u64 <= 2 * counts.w);
    let mut stats = cx.stats;
    stats.insert("groups", cx.groups.len() as u64);
    stats.insert("layers", nl as u64);
    stats.insert("added_cost", cx.cost as u64);
    Ok(Output { superstring: t, w: counts.w, stats })
}
