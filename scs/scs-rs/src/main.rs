//! scs-rs: polynomial-time 2-approximation for shortest common superstring
//! (openai/math #128) with compressed counts and lazily enabled periodic rules.
//!
//!   scs-rs <file> [-o out.txt] [--check]   one string per line, or FASTA
//!   scs-rs --batch <file>                  instances separated by '---' lines

mod approx;
mod counts;
mod sam;
mod verify;

use std::io::{BufRead, Write};
use std::time::Instant;

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
fn covers(t: &[u8], strings: &[Vec<u8>]) -> bool {
    let st = sam::Sam::build(vec![t.to_vec()]);
    strings.iter().all(|s| s.is_empty() || st.walk(s.iter().copied()).is_some())
}

fn read_strings(lines: &[String]) -> Vec<Vec<u8>> {
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let check = args.iter().any(|a| a == "--check");
    let batch = args.iter().any(|a| a == "--batch");
    let verify_lb = args.iter().any(|a| a == "--verify");
    let path = args.iter().skip(1).find(|a| !a.starts_with("--") && !a.starts_with('-')).expect("input file");
    let lines: Vec<String> = std::io::BufReader::new(std::fs::File::open(path).unwrap()).lines().map(|l| l.unwrap()).collect();
    if batch {
        let out = std::io::stdout();
        let mut out = out.lock();
        for chunk in lines.split(|l| l.trim() == "---") {
            let strings = read_strings(chunk);
            let r = solve(strings.clone(), check, verify_lb);
            assert!(covers(&r.superstring, &strings), "output misses an input");
            let st: Vec<String> = r.stats.iter().filter(|(k, v)| **v > 0 || k.starts_with("verify_")).map(|(k, v)| format!("{k}={v}")).collect();
            writeln!(out, "{} {} {} {} {}", r.w_ext, r.w, r.superstring.len(), r.rounds, st.join(",")).unwrap();
            writeln!(out, "{}", String::from_utf8_lossy(&r.superstring)).unwrap();
        }
        return;
    }
    let strings = read_strings(&lines);
    let total: usize = strings.iter().map(|s| s.len()).sum();
    let t0 = Instant::now();
    let r = solve(strings.clone(), check, verify_lb);
    let secs = t0.elapsed().as_secs_f64();
    assert!(covers(&r.superstring, &strings), "output misses an input");
    if let Some(i) = args.iter().position(|a| a == "-o") {
        std::fs::write(&args[i + 1], &r.superstring).unwrap();
    }
    let stats: Vec<String> = r.stats.iter().map(|(k, v)| format!("\"{k}\": {v}")).collect();
    println!(
        "{{\"strings\": {}, \"input_chars\": {}, \"seconds\": {:.2}, \"W_ext\": {}, \"W\": {}, \"length\": {}, \"length_over_W\": {:.4}, \"rounds\": {}, {}}}",
        strings.len(), total, secs, r.w_ext, r.w, r.superstring.len(),
        r.superstring.len() as f64 / r.w.max(1) as f64, r.rounds, stats.join(", ")
    );
}
