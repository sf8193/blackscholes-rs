//! Line-for-line port of Go `github.com/sf8193/blackscholes` v1.1.0
//! (`GetImpliedVolatility`, `BlackScholes`, `GetDelta`, `StdNormCDF`, `GetW`, `StdNormDensity`,
//! `GetGamma`, `GetVega`, `GetTheta`, `GetRho`).
//!
//! Parity matters: a Rust service replacing a Go one must return the same numbers. Price,
//! implied vol and delta match Go after its 4-decimal rounding; gamma, vega, theta and rho agree
//! within a relative 1e-12 (Go on amd64 uses an assembly `exp`). `libm` supplies `erfc`, `log`
//! and `exp`, and the arithmetic order is kept exactly as in Go.

/// `StdNormCDF`.
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * libm::erfc(-x / std::f64::consts::SQRT_2)
}

/// `GetW`.
pub fn w(s: f64, k: f64, t: f64, v: f64, r: f64) -> f64 {
    (r * t + v * v * t / 2.0 - libm::log(k / s)) / (v * t.sqrt())
}

/// `StdNormDensity`.
pub fn norm_density(x: f64) -> f64 {
    libm::exp(-(x * x) / 2.0) / (2.0 * std::f64::consts::PI).sqrt()
}

/// `BlackScholes` price of a call (`put == false`) or put.
pub fn price(s: f64, k: f64, t: f64, v: f64, r: f64, put: bool) -> f64 {
    let w = w(s, k, t, v, r);
    if !put {
        s * norm_cdf(w) - k * libm::exp(-r * t) * norm_cdf(w - v * t.sqrt())
    } else {
        k * libm::exp(-r * t) * norm_cdf(v * t.sqrt() - w) - s * norm_cdf(-w)
    }
}

/// `GetImpliedVolatility`: bisection, 100 iterations max, stops when cents (truncated) match.
pub fn implied_vol(cost: f64, s: f64, k: f64, t: f64, r: f64, put: bool, estimate: f64) -> f64 {
    let mut estimate = if estimate == 0.0 { 0.1 } else { estimate };
    let mut low = 0.0;
    let mut high = f64::INFINITY;
    for _ in 0..100 {
        let actual = price(s, k, t, estimate, r, put);
        if (cost * 100.0) as i64 == (actual * 100.0) as i64 {
            break;
        } else if actual > cost {
            high = estimate;
            estimate = (estimate - low) / 2.0 + low;
        } else {
            low = estimate;
            estimate = (high - estimate) / 2.0 + estimate;
            if !estimate.is_finite() {
                estimate = low * 2.0;
            }
        }
    }
    estimate
}

/// `GetDelta` for "CALL" (`put == false`) and "PUT".
pub fn delta(s: f64, k: f64, t: f64, v: f64, r: f64, put: bool) -> f64 {
    let w = w(s, k, t, v, r);
    let call = if !w.is_finite() {
        if s > k { 1.0 } else { 0.0 }
    } else {
        norm_cdf(w)
    };
    if !put {
        return call;
    }
    let d = call - 1.0;
    if d == -1.0 && k == s { 0.0 } else { d }
}

/// `GetGamma` (call and put alike).
pub fn gamma(s: f64, k: f64, t: f64, v: f64, r: f64) -> f64 {
    let w = w(s, k, t, v, r);
    if w.is_finite() {
        norm_density(w) / (s * v * t.sqrt())
    } else {
        0.0
    }
}

/// `GetVega` (call and put alike), per 1 vol point.
pub fn vega(s: f64, k: f64, t: f64, v: f64, r: f64) -> f64 {
    let w = w(s, k, t, v, r);
    if w.is_finite() {
        s * t.sqrt() * norm_density(w) / 100.0
    } else {
        0.0
    }
}

/// `GetTheta`; `scale` 0 means 365.
pub fn theta(s: f64, k: f64, t: f64, v: f64, r: f64, put: bool, scale: u32) -> f64 {
    let scale = if scale == 0 { 365 } else { scale };
    let w = w(s, k, t, v, r);
    if !w.is_finite() {
        return 0.0;
    }
    let decay = -v * s * norm_density(w) / (2.0 * t.sqrt());
    let carry = if !put {
        -(k * r * libm::exp(-r * t) * norm_cdf(w - v * t.sqrt()))
    } else {
        k * r * libm::exp(-r * t) * norm_cdf(v * t.sqrt() - w)
    };
    (decay + carry) / f64::from(scale)
}

/// `GetRho`; `scale` 0 means 100. Zero only when `w` is NaN (Go checks IsNaN, not isFinite).
pub fn rho(s: f64, k: f64, t: f64, v: f64, r: f64, put: bool, scale: u32) -> f64 {
    let scale = if scale == 0 { 100 } else { scale };
    let w = w(s, k, t, v, r);
    if w.is_nan() {
        return 0.0;
    }
    let x = if !put {
        k * t * libm::exp(-r * t) * norm_cdf(w - v * t.sqrt())
    } else {
        -k * t * libm::exp(-r * t) * norm_cdf(v * t.sqrt() - w)
    };
    x / f64::from(scale)
}
