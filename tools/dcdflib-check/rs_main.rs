#[allow(dead_code)]
#[path = "../../src/cdflib.rs"]
mod cdflib;
use std::io::{BufRead, Write};
/// Format like C's printf("%a").
fn hexf(v: f64) -> String {
    if v == 0.0 { return if v.is_sign_negative() { "-0x0p+0".into() } else { "0x0p+0".into() }; }
    let bits = v.to_bits();
    let sign = if bits >> 63 == 1 { "-" } else { "" };
    let exp = ((bits >> 52) & 0x7ff) as i64;
    let mut mant = bits & ((1u64 << 52) - 1);
    let (lead, e) = if exp == 0 { (0, -1022) } else { (1, exp - 1023) };
    let mut digits = 13;
    while digits > 0 && mant & 0xf == 0 { mant >>= 4; digits -= 1; }
    if digits == 0 { format!("{sign}0x{lead}p{e:+}") } else { format!("{sign}0x{lead}.{mant:0w$x}p{e:+}", w = digits) }
}
fn main() {
    let stdin = std::io::stdin();
    let mut out = std::io::BufWriter::new(std::io::stdout());
    let t = std::time::Instant::now();
    let mut n = 0;
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let f: Vec<&str> = line.split_whitespace().collect();
        let (u, v, w): (f64, f64, f64) = (f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap());
        if f[0] == "b" {
            let (c, cc) = cdflib::cumbet(w, 1.0 - w, u, v);
            writeln!(out, "{} {}", hexf(c), hexf(cc)).unwrap();
        } else {
            match cdflib::pnbinom(u, v, w) { Some(p) => writeln!(out, "{}", hexf(p)).unwrap(), None => writeln!(out, "undef").unwrap() }
        }
        n += 1;
    }
    eprintln!("{n} calls, {:.0} ns/call", t.elapsed().as_nanos() as f64 / n as f64);
}
