//! scs-rs: polynomial-time 2-approximation for shortest common superstring
//! (openai/math #128) with compressed counts and lazily enabled periodic rules.
//!
//!   scs-rs <file> [-o out.txt] [--check]   one string per line, or FASTA
//!   scs-rs --batch <file>                  instances separated by '---' lines

use lazy_scs_core::{covers, read_strings, solve};
use std::io::{BufRead, Write};
use std::time::Instant;

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
