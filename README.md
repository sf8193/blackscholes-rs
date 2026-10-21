# blackscholes

Black-Scholes price, implied volatility and greeks in Rust. It is a line-for-line port of
Go `github.com/sf8193/blackscholes` v1.1.0. Price, implied vol and delta match Go v1.1.0 after Go's 4-decimal rounding; gamma, vega, theta and
rho agree within a relative 1e-12 (Go on amd64 uses an assembly exp; this crate uses libm).

Functions: `norm_cdf`, `norm_density`, `w`, `price`, `implied_vol`, `delta`, `gamma`, `vega`,
`theta`, `rho`. Only dependency: `libm`.

Go v1.1.0 `GetTheta`/`GetRho` test lowercase `"call"` (`GetDelta` tests `"CALL"`); here `put: bool`
means call or put as it should, and the oracle passes `"call"`/`"put"` to match.

## Golden test

`tests/fixtures/bs_golden.json` was made by running the Go library itself over a fixed
random set of inputs (`oracle/main.go`, seed 42); `tests/fixtures/greeks_golden.json` is the same for
gamma/vega/theta/rho (seed 43). `cargo test` checks every row against them.

Regenerate on linux/amd64:

    docker run --rm --platform linux/amd64 -v "$PWD/oracle:/w" -w /w golang:1.23 \
      sh -c 'go run . out.json greeks.json'

then copy `oracle/out.json` to `tests/fixtures/bs_golden.json` and `oracle/greeks.json` to
`tests/fixtures/greeks_golden.json`.
