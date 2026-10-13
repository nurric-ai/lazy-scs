//! Independent certificate check for the lower bound W' <= OPT.
//!
//! Every count m'(s) must be justified by bounds that only mention strictly longer
//! words: an extension sum (2.5), the required-string bound, or a periodic rule (2.6)
//! on an enabled text whose trigger holds. If every count is justified, induction on
//! decreasing length (the paper's Lemma 2.1, generalised to any rule subset) gives
//! m' <= N_T for every common superstring T, hence W' <= OPT.
//!
//! This module recomputes all of that from the automaton and the final step
//! functions only. It does not reuse the solver's stored R tables, children sums,
//! rule list or blocking-functional code.

use crate::counts::Counts;
use crate::sam::{key, FxMap, Sam};

fn m_at(c: &Counts, x: u32, l: u32) -> u64 {
    let f = c.f(x);
    let mut v = 0;
    for &(st, val) in f {
        if st <= l {
            v = val;
        }
    }
    v as u64
}

fn m_word(sam: &Sam, c: &Counts, w: &[u8]) -> u64 {
    match sam.walk(w.iter().copied()) {
        Some(x) if !w.is_empty() => m_at(c, x, w.len() as u32),
        _ => 0,
    }
}

/// Straightforward blocking functional (2.3): A = rotation `r` of `g` repeated,
/// s = A[0:slen]. Every membership test walks from the root.
fn blocking(sam: &Sam, c: &Counts, g: &[u8], r: usize, slen: usize, lm: usize, sigma: &[u8]) -> u64 {
    let p = g.len() as i64;
    let a_at = |i: i64| g[((i + r as i64).rem_euclid(p)) as usize];
    let word = |lo: i64, hi: i64| -> Vec<u8> { (lo..hi).map(a_at).collect() };
    let mut total = 0u64;
    let mut a = 0i64;
    while a as usize + slen <= lm && sam.walk(word(-a, slen as i64)).is_some() {
        let mut cc = 0i64;
        while a as usize + slen + cc as usize <= lm {
            let m = word(-a, slen as i64 + cc);
            if sam.walk(m.iter().copied()).is_none() {
                break;
            }
            let (xbad, ybad) = (a_at(-a - 1), a_at(slen as i64 + cc));
            for &x in sigma {
                for &y in sigma {
                    if x != xbad && y != ybad {
                        let mut r = vec![x];
                        r.extend_from_slice(&m);
                        r.push(y);
                        total += m_word(sam, c, &r);
                    }
                }
            }
            cc += 1;
        }
        a += 1;
    }
    total
}

pub struct Report {
    pub checked_points: u64,
    pub unjustified: u64,
    pub ext_violations: u64,
}

pub fn verify(sam: &Sam, c: &Counts, texts: &[(usize, Vec<u8>)]) -> Report {
    let n = sam.n();
    let lm = sam.strings.iter().map(|s| s.len()).max().unwrap_or(0);
    let mut sigma: Vec<u8> = sam.next[0].iter().map(|e| e.0).collect();
    sigma.sort_unstable();

    // children sums recomputed from the suffix-link tree
    let mut child = vec![0u64; n];
    for z in 1..n as u32 {
        let lk = sam.link[z as usize];
        let minlen = sam.len[lk as usize] + 1;
        child[lk as usize] += m_at(c, z, minlen);
    }
    let mut is_req = vec![false; n];
    for s in &sam.strings {
        is_req[sam.walk(s.iter().copied()).unwrap() as usize] = true;
    }

    // periodic rule bounds, recomputed: (state, len) of w -> best triggered bound
    let mut rule_bound: FxMap<u64, u64> = FxMap::default();
    for (p, g) in texts {
        let p = *p;
        for r in 0..p {
            let mut path = vec![0u32];
            let mut v = 0u32;
            let mut i = 0;
            while let Some(nv) = sam.get(v, g[(r + i) % p]) {
                v = nv;
                path.push(v);
                i += 1;
            }
            for nlen in p + 1..path.len() {
                let trig = m_at(c, path[nlen], nlen as u32) > blocking(sam, c, g, r, nlen, lm, &sigma);
                if !trig {
                    continue;
                }
                for k in 1..=(nlen - 1) / p {
                    let wl = nlen - k * p;
                    let b = blocking(sam, c, g, r, wl, lm, &sigma) + k as u64 + 1;
                    let e = rule_bound.entry(key(path[wl], wl as u32)).or_insert(0);
                    *e = (*e).max(b);
                }
            }
        }
    }

    let mut rep = Report { checked_points: 0, unjustified: 0, ext_violations: 0 };
    for x in 1..n as u32 {
        let xi = x as usize;
        let lo = sam.len[sam.link[xi] as usize] + 1;
        let hi = sam.len[xi];
        // only the last length of each constant piece of m can need more than m(l+1)
        let mut pts: Vec<u32> = c.f(x).iter().filter(|s| s.0 > lo).map(|s| s.0 - 1).collect();
        pts.push(hi);
        pts.extend(c.bumps.keys().filter(|&&k| (k >> 32) as u32 == x).map(|&k| k as u32));
        // R changes inside a piece too: check where any child's function steps
        for &(_, y) in &sam.next[xi] {
            pts.extend(c.f(y).iter().filter(|s| s.0 >= 1 && s.0 - 1 >= lo && s.0 - 1 <= hi).map(|s| s.0 - 1));
        }
        pts.sort_unstable();
        pts.dedup();
        for l in pts {
            rep.checked_points += 1;
            let m = m_at(c, x, l);
            let left = if l < hi { m_at(c, x, l + 1) } else { child[xi] };
            let right: u64 = sam.next[xi].iter().map(|&(_, y)| m_at(c, y, l + 1)).sum();
            if m < left || m < right {
                rep.ext_violations += 1;
            }
            let req = if l == hi && is_req[xi] { 1 } else { 0 };
            let rule = *rule_bound.get(&key(x, l)).unwrap_or(&0);
            if m > left.max(right).max(req).max(rule) {
                rep.unjustified += 1;
            }
        }
    }
    rep
}
