//! Forced occurrence counts m(s) as per-state step functions (port of compressed.py).
//! Extension rules (2.5) + required strings exactly; periodic rules (2.6) only for
//! enabled texts, applied as point lower bounds iterated to a fixpoint.

use crate::sam::{key, kl, kx, FxMap, Sam, NONE};

#[inline]
pub fn ev(segs: &[(u32, u32)], l: u32) -> u32 {
    let i = segs.partition_point(|s| s.0 <= l);
    segs[i - 1].1
}

pub struct Counts {
    f_segs: Vec<(u32, u32)>,
    f_idx: Vec<(u32, u32)>,
    r_segs: Vec<(u32, u32)>,
    r_idx: Vec<(u32, u32)>,
    pub children_sum: Vec<u32>,
    pub is_req: Vec<bool>,
    pub w: u64,
    pub bumps: FxMap<u64, u32>,
    pub passes: usize,
}

impl Counts {
    pub fn f(&self, x: u32) -> &[(u32, u32)] {
        let (o, c) = self.f_idx[x as usize];
        &self.f_segs[o as usize..(o + c) as usize]
    }
    pub fn r(&self, x: u32) -> &[(u32, u32)] {
        let (o, c) = self.r_idx[x as usize];
        &self.r_segs[o as usize..(o + c) as usize]
    }
    #[inline]
    pub fn m(&self, k: u64) -> u32 {
        ev(self.f(kx(k)), kl(k))
    }

    pub fn new(sam: &Sam, order: &[u32], texts: &[(usize, Vec<u8>)]) -> Counts {
        let n = sam.n();
        let mut is_req = vec![false; n];
        for s in &sam.strings {
            is_req[sam.walk(s.iter().copied()).unwrap() as usize] = true;
        }
        let mut c = Counts {
            f_segs: Vec::new(),
            f_idx: vec![(0, 0); n],
            r_segs: Vec::new(),
            r_idx: vec![(0, 0); n],
            children_sum: vec![0; n],
            is_req,
            w: 0,
            bumps: FxMap::default(),
            passes: 0,
        };
        c.compute(sam, order);
        if texts.is_empty() {
            return c;
        }
        let maxlen = sam.strings.iter().map(|s| s.len()).max().unwrap_or(0);
        let mut sigma: Vec<u8> = sam.next[0].iter().map(|e| e.0).collect();
        sigma.sort_unstable();
        let rules = rules(sam, texts);
        for pass in 1..=100_000 {
            c.passes = pass;
            let mut changed = false;
            let mut cache: FxMap<(usize, usize, u32), u64> = FxMap::default();
            for r in &rules {
                let tx = &texts[r.ti].1;
                let bv = c.blocking(sam, tx, r.ti, r.pi, r.vl, maxlen, &sigma, &mut cache);
                if c.m(key(r.vx, r.vl)) as u64 > bv {
                    let need = c.blocking(sam, tx, r.ti, r.pi, r.wl, maxlen, &sigma, &mut cache) + r.k + 1;
                    if (c.m(key(r.wx, r.wl)) as u64) < need {
                        c.bumps.insert(key(r.wx, r.wl), need as u32);
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
            c.compute(sam, order);
        }
        c
    }

    fn compute(&mut self, sam: &Sam, order: &[u32]) {
        let n = sam.n();
        self.f_segs.clear();
        self.r_segs.clear();
        self.children_sum = vec![0; n];
        let mut bump_at: FxMap<u32, Vec<u32>> = FxMap::default();
        for &k in self.bumps.keys() {
            bump_at.entry(kx(k)).or_default().push(kl(k));
        }
        let mut pts: Vec<u32> = Vec::new();
        let mut r: Vec<(u32, u32)> = Vec::new();
        let mut segs: Vec<(u32, u32)> = Vec::new();
        for &x in order {
            if x == 0 {
                continue;
            }
            let xi = x as usize;
            let lo = sam.len[sam.link[xi] as usize] + 1;
            let hi = sam.len[xi];
            pts.clear();
            pts.push(lo);
            for &(_, y) in &sam.next[xi] {
                for &(st, _) in self.f(y) {
                    if st >= 1 && st - 1 >= lo && st - 1 <= hi {
                        pts.push(st - 1);
                    }
                }
            }
            pts.sort_unstable();
            pts.dedup();
            r.clear();
            for &l in &pts {
                let mut s = 0u32;
                for &(_, y) in &sam.next[xi] {
                    s += ev(self.f(y), l + 1);
                }
                r.push((l, s));
            }
            let ro = self.r_segs.len() as u32;
            self.r_segs.extend_from_slice(&r);
            self.r_idx[xi] = (ro, r.len() as u32);

            let mut evalpts: Vec<(u32, u32)>;
            let bl = bump_at.get(&x);
            if let Some(bl) = bl {
                let mut cut: Vec<u32> = pts.clone();
                for &l in bl {
                    cut.push(l);
                    if l < hi {
                        cut.push(l + 1);
                    }
                }
                cut.sort_unstable();
                cut.dedup();
                evalpts = cut.iter().map(|&l| (l, ev(&r, l))).collect();
            } else {
                evalpts = std::mem::take(&mut r);
            }
            let mut cur = self.children_sum[xi].max(self.is_req[xi] as u32);
            segs.clear();
            for &(l, val) in evalpts.iter().rev() {
                let v = if bl.is_some() { val.max(*self.bumps.get(&key(x, l)).unwrap_or(&0)) } else { val };
                cur = cur.max(v);
                segs.push((l, cur));
            }
            segs.reverse();
            let fo = self.f_segs.len() as u32;
            for &(l, v) in &segs {
                if self.f_segs.len() as u32 == fo || self.f_segs.last().unwrap().1 != v {
                    self.f_segs.push((l, v));
                }
            }
            self.f_idx[xi] = (fo, self.f_segs.len() as u32 - fo);
            let first = self.f_segs[fo as usize].1;
            self.children_sum[sam.link[xi] as usize] += first;
            if bl.is_none() {
                r = std::mem::take(&mut evalpts);
            }
        }
        self.w = sam.next[0].iter().map(|&(_, y)| ev(self.f(y), 1) as u64).sum();
    }

    /// Blocking functional (2.3) for s = A[0:slen], A = rotation pi of P repeated.
    fn blocking(&self, sam: &Sam, p: &[u8], ti: usize, pi: usize, slen: u32, lm: usize, sigma: &[u8],
                cache: &mut FxMap<(usize, usize, u32), u64>) -> u64 {
        let ck = (ti, pi, slen);
        if let Some(&v) = cache.get(&ck) {
            return v;
        }
        let pl = p.len();
        let slen = slen as usize;
        let off = lm + 1;
        let span = off + slen + lm + 2;
        let st = (pi + pl - (off % pl)) % pl;
        let a_str: Vec<u8> = (0..span).map(|i| p[(st + i) % pl]).collect();
        // M = A[-a : slen + c]; states of M and of xM are extended one letter at a time
        let mut total = 0u64;
        let mut a = 0usize;
        while a + slen <= lm {
            let m0 = match sam.walk(a_str[off - a..off + slen].iter().copied()) {
                Some(v) => v,
                None => break,
            };
            let xbad = a_str[off - a - 1];
            let mut xs: Vec<u32> = sigma
                .iter()
                .filter(|&&x| x != xbad)
                .filter_map(|&x| sam.get(0, x).and_then(|s0| {
                    let mut v = s0;
                    for &ch in &a_str[off - a..off + slen] {
                        v = sam.get(v, ch)?;
                    }
                    Some(v)
                }))
                .collect();
            let mut vm = m0;
            let mut mlen = a + slen;
            let mut c = 0usize;
            loop {
                let ybad = a_str[off + slen + c];
                for &xm in &xs {
                    for &y in sigma {
                        if y != ybad {
                            if let Some(z) = sam.get(xm, y) {
                                total += ev(self.f(z), (mlen + 2) as u32) as u64;
                            }
                        }
                    }
                }
                if a + slen + c + 1 > lm {
                    break;
                }
                let ch = a_str[off + slen + c];
                vm = match sam.get(vm, ch) {
                    Some(v) => v,
                    None => break,
                };
                xs = xs.iter().filter_map(|&s| sam.get(s, ch)).collect();
                c += 1;
                mlen += 1;
            }
            let _ = vm;
            a += 1;
        }
        cache.insert(ck, total);
        total
    }

    /// Nonzero edge multiplicities: up[w] into word w, down[w] out of w.
    pub fn edges(&self, sam: &Sam) -> (FxMap<u64, u32>, FxMap<u64, u32>) {
        let mut up: FxMap<u64, u32> = FxMap::default();
        let mut down: FxMap<u64, u32> = FxMap::default();
        for x in 1..sam.n() as u32 {
            let xi = x as usize;
            let lo = sam.len[sam.link[xi] as usize] + 1;
            let hi = sam.len[xi];
            let fx = self.f(x);
            let rx = self.r(x);
            let u_hi = ev(fx, hi) as i64 - self.children_sum[xi] as i64;
            assert!(u_hi >= 0);
            if u_hi > 0 {
                up.insert(key(x, hi), u_hi as u32);
            }
            for &(st, val) in &fx[1..] {
                let u = ev(fx, st - 1) - val;
                if u > 0 {
                    up.insert(key(x, st - 1), u);
                }
            }
            let mut cuts: Vec<u32> = fx.iter().map(|s| s.0).chain(rx.iter().map(|s| s.0)).collect();
            cuts.sort_unstable();
            cuts.dedup();
            for (i, &a) in cuts.iter().enumerate() {
                let b = if i + 1 < cuts.len() { cuts[i + 1] - 1 } else { hi };
                let fv = ev(fx, a);
                let rv = ev(rx, a);
                assert!(fv >= rv);
                if fv > rv {
                    for l in a..=b {
                        down.insert(key(x, l), fv - rv);
                    }
                }
            }
            let _ = lo;
        }
        let _ = NONE;
        (up, down)
    }
}

pub struct Rule {
    vx: u32,
    vl: u32,
    wx: u32,
    wl: u32,
    ti: usize,
    pi: usize,
    k: u64,
}

/// All rules (v, p, k) with v on an enabled periodic text, as (state, len) pairs.
fn rules(sam: &Sam, texts: &[(usize, Vec<u8>)]) -> Vec<Rule> {
    let mut out = Vec::new();
    for (ti, (p, g)) in texts.iter().enumerate() {
        let p = *p;
        for r in 0..p {
            let mut path = vec![0u32];
            let mut v = 0u32;
            let mut i = 0usize;
            while let Some(nv) = sam.get(v, g[(r + i) % p]) {
                v = nv;
                path.push(v);
                i += 1;
            }
            for n in p + 1..path.len() {
                for k in 1..=(n - 1) / p {
                    out.push(Rule {
                        vx: path[n],
                        vl: n as u32,
                        wx: path[n - k * p],
                        wl: (n - k * p) as u32,
                        ti,
                        pi: r,
                        k: k as u64,
                    });
                }
            }
        }
    }
    out
}
