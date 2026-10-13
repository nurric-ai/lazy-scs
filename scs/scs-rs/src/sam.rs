//! Generalized suffix automaton. A word is identified by (state, length).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

pub const NONE: u32 = u32::MAX;

#[derive(Default)]
pub struct Fx(u64);
impl Hasher for Fx {
    fn write(&mut self, b: &[u8]) {
        for &x in b {
            self.0 = (self.0.rotate_left(5) ^ x as u64).wrapping_mul(0x517cc1b727220a95);
        }
    }
    fn write_u64(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x517cc1b727220a95);
    }
    fn write_usize(&mut self, i: usize) {
        self.write_u64(i as u64)
    }
    fn finish(&self) -> u64 {
        // splitmix64 finalizer: spread high-bit entropy into the low bits used for buckets
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
}
pub type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;

#[inline]
pub fn key(x: u32, l: u32) -> u64 {
    ((x as u64) << 32) | l as u64
}
#[inline]
pub fn kx(k: u64) -> u32 {
    (k >> 32) as u32
}
#[inline]
pub fn kl(k: u64) -> u32 {
    k as u32
}
pub const EPS: u64 = 0;

pub struct Sam {
    pub next: Vec<Vec<(u8, u32)>>,
    pub link: Vec<u32>,
    pub len: Vec<u32>,
    pub pos: Vec<(u32, u32)>, // one occurrence: (string index, end)
    pub cnt: Vec<u32>,        // |endpos| after propagate()
    pub strings: Vec<Vec<u8>>,
}

impl Sam {
    pub fn build(strings: Vec<Vec<u8>>) -> Sam {
        let mut s = Sam {
            next: vec![Vec::new()],
            link: vec![NONE],
            len: vec![0],
            pos: vec![(0, 0)],
            cnt: vec![0],
            strings: Vec::new(),
        };
        for (si, st) in strings.iter().enumerate() {
            let mut last = 0u32;
            for (i, &c) in st.iter().enumerate() {
                last = s.extend(last, c, (si as u32, i as u32 + 1));
                s.cnt[last as usize] += 1;
            }
        }
        s.strings = strings;
        s
    }

    pub fn n(&self) -> usize {
        self.len.len()
    }

    #[inline]
    pub fn get(&self, v: u32, c: u8) -> Option<u32> {
        self.next[v as usize].iter().find(|e| e.0 == c).map(|e| e.1)
    }

    fn set(&mut self, v: u32, c: u8, t: u32) {
        let row = &mut self.next[v as usize];
        if let Some(e) = row.iter_mut().find(|e| e.0 == c) {
            e.1 = t
        } else {
            row.push((c, t));
            row.sort_unstable();
        }
    }

    fn new_state(&mut self, len: u32, next: Vec<(u8, u32)>, link: u32, pos: (u32, u32)) -> u32 {
        self.next.push(next);
        self.link.push(link);
        self.len.push(len);
        self.pos.push(pos);
        self.cnt.push(0);
        (self.len.len() - 1) as u32
    }

    fn split(&mut self, p: u32, c: u8, q: u32) -> u32 {
        let qi = q as usize;
        let clone = self.new_state(self.len[p as usize] + 1, self.next[qi].clone(), self.link[qi], self.pos[qi]);
        let mut p = p;
        while p != NONE && self.get(p, c) == Some(q) {
            self.set(p, c, clone);
            p = self.link[p as usize];
        }
        self.link[qi] = clone;
        clone
    }

    fn extend(&mut self, last: u32, c: u8, pos: (u32, u32)) -> u32 {
        if let Some(q) = self.get(last, c) {
            return if self.len[last as usize] + 1 == self.len[q as usize] { q } else { self.split(last, c, q) };
        }
        let cur = self.new_state(self.len[last as usize] + 1, Vec::new(), NONE, pos);
        let mut p = last;
        while p != NONE && self.get(p, c).is_none() {
            self.set(p, c, cur);
            p = self.link[p as usize];
        }
        if p == NONE {
            self.link[cur as usize] = 0;
        } else {
            let q = self.get(p, c).unwrap();
            self.link[cur as usize] =
                if self.len[p as usize] + 1 == self.len[q as usize] { q } else { self.split(p, c, q) };
        }
        cur
    }

    /// States ordered by decreasing length (root last).
    pub fn order_desc(&self) -> Vec<u32> {
        let maxl = *self.len.iter().max().unwrap() as usize;
        let mut bucket = vec![0usize; maxl + 2];
        for &l in &self.len {
            bucket[l as usize] += 1;
        }
        let mut start = vec![0usize; maxl + 2];
        let mut acc = 0;
        for l in (0..=maxl).rev() {
            start[l] = acc;
            acc += bucket[l];
        }
        let mut out = vec![0u32; self.n()];
        for (x, &l) in self.len.iter().enumerate() {
            out[start[l as usize]] = x as u32;
            start[l as usize] += 1;
        }
        out
    }

    /// cnt[x] = |endpos(x)|.
    pub fn propagate(&mut self) {
        for x in self.order_desc() {
            let l = self.link[x as usize];
            if l != NONE {
                self.cnt[l as usize] += self.cnt[x as usize];
            }
        }
    }

    pub fn walk<I: IntoIterator<Item = u8>>(&self, s: I) -> Option<u32> {
        let mut v = 0u32;
        for c in s {
            v = self.get(v, c)?;
        }
        Some(v)
    }

    /// Vertex of a word given as bytes.
    pub fn vertex<I: IntoIterator<Item = u8>>(&self, s: I) -> Option<u64> {
        let mut v = 0u32;
        let mut l = 0u32;
        for c in s {
            v = self.get(v, c)?;
            l += 1;
        }
        Some(key(v, l))
    }

    #[inline]
    pub fn suf(&self, k: u64) -> u64 {
        let (x, l) = (kx(k), kl(k));
        debug_assert!(l > 0);
        let lk = self.link[x as usize];
        if l - 1 > self.len[lk as usize] { key(x, l - 1) } else { key(lk, l - 1) }
    }

    #[inline]
    pub fn ext(&self, k: u64, c: u8) -> Option<u64> {
        self.get(kx(k), c).map(|y| key(y, kl(k) + 1))
    }

    #[inline]
    pub fn last_char(&self, x: u32) -> u8 {
        let (si, e) = self.pos[x as usize];
        self.strings[si as usize][e as usize - 1]
    }

    pub fn word(&self, k: u64) -> &[u8] {
        let (si, e) = self.pos[kx(k) as usize];
        &self.strings[si as usize][(e - kl(k)) as usize..e as usize]
    }
}

/// Single-string automaton plus, per state, the largest end position of its words.
pub fn with_max_end(t: Vec<u8>) -> (Sam, Vec<u32>) {
    let mut s = Sam { next: vec![Vec::new()], link: vec![NONE], len: vec![0], pos: vec![(0, 0)], cnt: vec![0], strings: Vec::new() };
    let mut lasts = Vec::with_capacity(t.len());
    let mut last = 0u32;
    for (i, &c) in t.iter().enumerate() {
        last = s.extend(last, c, (0, i as u32 + 1));
        lasts.push(last);
    }
    let mut max_end = vec![0u32; s.n()];
    for (i, &x) in lasts.iter().enumerate() {
        max_end[x as usize] = i as u32 + 1;
    }
    for x in s.order_desc() {
        let l = s.link[x as usize];
        if l != NONE {
            max_end[l as usize] = max_end[l as usize].max(max_end[x as usize]);
        }
    }
    s.strings = vec![t];
    (s, max_end)
}

/// Drop leading then trailing letters while every input still occurs (never lengthens).
pub fn trim(t: &[u8], req: &[Vec<u8>]) -> Vec<u8> {
    if req.is_empty() {
        return Vec::new();
    }
    let (s1, max_end) = with_max_end(t.to_vec());
    let k = req.iter().map(|w| max_end[s1.walk(w.iter().copied()).unwrap() as usize] as usize - w.len()).min().unwrap();
    let t2 = &t[k..];
    let s2 = Sam::build(vec![t2.to_vec()]);
    let e = req.iter().map(|w| s2.pos[s2.walk(w.iter().copied()).unwrap() as usize].1 as usize).max().unwrap();
    t2[..e].to_vec()
}

/// Drop duplicates and strings contained in another (OPT-preserving); sorted.
pub fn reduce_inputs(strings: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut u: Vec<Vec<u8>> = strings.into_iter().filter(|s| !s.is_empty()).collect();
    u.sort_unstable();
    u.dedup();
    let mut sam = Sam::build(u);
    sam.propagate();
    let mut keep: Vec<Vec<u8>> = Vec::new();
    for s in &sam.strings {
        let x = sam.walk(s.iter().copied()).unwrap();
        if sam.cnt[x as usize] == 1 {
            keep.push(s.clone());
        }
    }
    keep
}
