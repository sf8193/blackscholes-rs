use blackscholes::{delta, implied_vol};
use serde::Deserialize;

#[derive(Deserialize)]
struct Row {
    #[serde(rename = "S")]
    s: f64,
    #[serde(rename = "K")]
    k: f64,
    #[serde(rename = "T")]
    t: f64,
    #[serde(rename = "Mark")]
    mark: f64,
    #[serde(rename = "R")]
    r: f64,
    #[serde(rename = "Put")]
    put: bool,
    iv_bits: u64,
    delta_bits: u64,
    iv_r4: f64,
    delta_r4: f64,
    est: f64,
    v: Option<f64>,
}

fn round4(x: f64) -> f64 {
    (x * 10000.0).round() / 10000.0
}

/// Rust against Go's own output (sf8193/blackscholes v1.1.0), bit for bit.
/// Random rows may differ in the last bits if round4 agrees; the estimate-0 and fixed-vol rows
/// must be bit-exact.
#[test]
fn matches_go_oracle() {
    let rows: Vec<Row> = serde_json::from_str(include_str!("fixtures/bs_golden.json")).unwrap();
    let (mut exact, mut r4, mut bad) = (0, 0, vec![]);
    for (i, w) in rows.iter().enumerate() {
        let iv =
            w.v.unwrap_or_else(|| implied_vol(w.mark, w.s, w.k, w.t, w.r, w.put, w.est));
        let d = delta(w.s, w.k, w.t, iv, w.r, w.put);
        if iv.to_bits() == w.iv_bits && d.to_bits() == w.delta_bits {
            exact += 1;
        } else if w.est == 0.11
            && w.v.is_none()
            && round4(iv * 100.0) == w.iv_r4
            && round4(d) == w.delta_r4
        {
            r4 += 1;
        } else {
            bad.push(i);
        }
    }
    println!(
        "rows {} bit-exact {exact} round4-only {r4} differ {}",
        rows.len(),
        bad.len()
    );
    assert!(bad.is_empty(), "differ after round4: {bad:?}");
    assert_eq!(rows.len(), exact + r4);
    assert_eq!(
        r4, 211,
        "bit-exact count changed: libm or the arithmetic order moved"
    );
}

#[derive(Deserialize)]
struct GRow {
    #[serde(rename = "S")]
    s: f64,
    #[serde(rename = "K")]
    k: f64,
    #[serde(rename = "T")]
    t: f64,
    #[serde(rename = "V")]
    v: f64,
    #[serde(rename = "R")]
    r: f64,
    #[serde(rename = "Put")]
    put: bool,
    #[serde(rename = "Scale")]
    scale: u32,
    #[serde(rename = "GammaBits")]
    gamma: u64,
    #[serde(rename = "VegaBits")]
    vega: u64,
    #[serde(rename = "ThetaBits")]
    theta: u64,
    #[serde(rename = "RhoBits")]
    rho: u64,
}

/// gamma, vega, theta and rho against Go v1.1.0: relative difference <= 1e-12 on every row.
/// Not bit-exact: Go on amd64 uses an assembly `exp`, this crate uses `libm`.
#[test]
fn greeks_match_go_oracle() {
    let rows: Vec<GRow> =
        serde_json::from_str(include_str!("fixtures/greeks_golden.json")).unwrap();
    let (mut exact, mut bad) = ([0usize; 4], [0usize; 4]);
    for g in &rows {
        let got = [
            blackscholes::gamma(g.s, g.k, g.t, g.v, g.r),
            blackscholes::vega(g.s, g.k, g.t, g.v, g.r),
            blackscholes::theta(g.s, g.k, g.t, g.v, g.r, g.put, g.scale),
            blackscholes::rho(g.s, g.k, g.t, g.v, g.r, g.put, g.scale),
        ];
        for (i, (x, want)) in got
            .iter()
            .zip([g.gamma, g.vega, g.theta, g.rho])
            .enumerate()
        {
            let w = f64::from_bits(want);
            if x.to_bits() == want {
                exact[i] += 1;
            }
            let rel = (x - w).abs() / w.abs().max(1e-300);
            if rel.is_nan() || rel > 1e-12 {
                bad[i] += 1;
            }
        }
    }
    println!(
        "rows {} bit-exact [gamma, vega, theta, rho] = {exact:?} over 1e-12 = {bad:?}",
        rows.len()
    );
    assert_eq!(bad, [0; 4]);
}
