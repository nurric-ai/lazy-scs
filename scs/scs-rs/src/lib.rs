//! lazy-scs core: polynomial-time 2-approximation for shortest common superstring
//! (openai/math #128) with compressed counts and lazily enabled periodic rules.
//! Used by the `scs-rs` command line (src/main.rs) and, with the `python` feature, by the
//! `lazy_scs` Python package (src/python.rs).

mod approx;
mod counts;
mod sam;
mod verify;
#[cfg(feature = "python")]
mod python;

pub struct Solved {
    pub superstring: Vec<u8>,
    pub w_ext: u64,
    pub w: u64,
    pub rounds: usize,
    pub stats: std::collections::BTreeMap<&'static str, u64>,
}

pub fn solve(strings: Vec<Vec<u8>>, check: bool, verify_lb: bool) -> Solved {
    let req = sam::reduce_inputs(strings);
    if req.is_empty() {
        return Solved { superstring: vec![], w_ext: 0, w: 0, rounds: 0, stats: Default::default() };
    }
    let s = sam::Sam::build(req);
    let order = s.order_desc();
    let mut texts: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut w_ext = 0;
    for round in 1.. {
        let c = counts::Counts::new(&s, &order, &texts);
        if round == 1 {
            w_ext = c.w;
        }
        match approx::run(&s, &c, check) {
            Ok(o) => {
                let mut stats = o.stats;
                stats.insert("count_passes", c.passes as u64);
                if verify_lb {
                    let rep = verify::verify(&s, &c, &texts);
                    stats.insert("verify_points", rep.checked_points);
                    stats.insert("verify_unjustified", rep.unjustified);
                    stats.insert("verify_ext_violations", rep.ext_violations);
                }
                stats.insert("untrimmed_length", o.superstring.len() as u64);
                let trimmed = sam::trim(&o.superstring, &s.strings);
                return Solved { superstring: trimmed, w_ext, w: o.w, rounds: round, stats };
            }
            Err(approx::NeedRule(p, g)) => {
                assert!(!texts.iter().any(|t| t.0 == p && t.1 == g), "rule enabled but still missing");
                texts.push((p, g));
            }
        }
    }
    unreachable!()
}

/// Keep every input: verify with a suffix automaton of the output.
pub fn covers(t: &[u8], strings: &[Vec<u8>]) -> bool {
    let st = sam::Sam::build(vec![t.to_vec()]);
    strings.iter().all(|s| s.is_empty() || st.walk(s.iter().copied()).is_some())
}

/// One string per line, or FASTA records if any line starts with '>'.
pub fn read_strings(lines: &[String]) -> Vec<Vec<u8>> {
    if lines.iter().any(|l| l.starts_with('>')) {
        let mut out = Vec::new();
        let mut cur = Vec::new();
        for l in lines {
            if l.starts_with('>') {
                if !cur.is_empty() {
                    out.push(std::mem::take(&mut cur));
                }
            } else {
                cur.extend(l.trim().bytes());
            }
        }
        if !cur.is_empty() {
            out.push(cur);
        }
        out
    } else {
        lines.iter().map(|l| l.trim_end_matches(['\r', '\n']).as_bytes().to_vec()).filter(|s| !s.is_empty()).collect()
    }
}

