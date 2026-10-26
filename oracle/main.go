package main

import (
	"encoding/json"
	"math"
	"math/rand"
	"os"

	"github.com/sf8193/blackscholes"
)

type Row struct {
	S, K, T, Mark, R float64
	Put              bool
	IVBits           uint64   `json:"iv_bits"`
	DeltaBits        uint64   `json:"delta_bits"`
	IVR4             float64  `json:"iv_r4"`
	DeltaR4          float64  `json:"delta_r4"`
	Est              float64  `json:"est"`
	V                *float64 `json:"v"` // fixed vol: skip the solver, delta only
}

func round4(v float64) float64 { return math.Round(v*10000) / 10000 }

func main() {
	rng := rand.New(rand.NewSource(42))
	const r = 0.0425
	var rows []Row
	addE := func(s, k, t, mark float64, put bool, est float64, fixedV *float64) {
		ct := "CALL"
		if put {
			ct = "PUT"
		}
		iv := blackscholes.GetImpliedVolatility(mark, s, k, t, r, ct, est)
		if fixedV != nil {
			iv = *fixedV
		}
		d, err := blackscholes.GetDelta(s, k, t, iv, r, ct)
		if err != nil {
			panic(err)
		}
		rows = append(rows, Row{s, k, t, mark, r, put, math.Float64bits(iv), math.Float64bits(d), round4(iv * 100), round4(d), est, fixedV})
	}
	add := func(s, k, t, mark float64, put bool) { addE(s, k, t, mark, put, 0.11, nil) }
	for i := 0; i < 4000; i++ {
		s := 7700.0
		if i%2 == 1 {
			s = 770
		}
		s += (rng.Float64() - 0.5) * 0.02 * s
		k := math.Round(s*(1+(rng.Float64()*0.3-0.15))/5) * 5
		t := math.Exp(math.Log(1.0/(365.25*24*60)) + rng.Float64()*(math.Log(21.0/365.25)-math.Log(1.0/(365.25*24*60))))
		mark := math.Exp(math.Log(0.05) + rng.Float64()*(math.Log(400)-math.Log(0.05)))
		mark = math.Round(mark*100) / 100
		if mark < 0.05 {
			mark = 0.05
		}
		put := i%4 >= 2
		if i%3 != 0 { // realistic: mark from Go's own price at a random vol, in cents
			ct := "CALL"
			if put {
				ct = "PUT"
			}
			mark = math.Round(blackscholes.BlackScholes(s, k, t, 0.08+rng.Float64()*0.5, r, ct)*100) / 100
			if mark < 0.05 {
				mark = 0.05
			}
		}
		add(s, k, t, mark, put)
	}
	// edge rows
	for _, put := range []bool{false, true} {
		for _, s := range []float64{7700, 770} {
			day := 1.0 / 365.25
			minute := day / 1440
			add(s, s, day, 5, put)          // k == s
			add(s, s, minute, 0.05, put)    // tiny T, at the money
			add(s, s, minute/60, 0.05, put) // one second
			add(s, s*0.7, day, 0.05, put)   // deep OTM/ITM
			add(s, s*1.3, day, 0.05, put)
			add(s, s*0.7, day, 400, put)
			add(s, s*1.3, day, 400, put)
			add(s, s*0.9, minute, 0.05, put)
			add(s, s*1.1, minute, 0.05, put)
			add(s, s*0.95, day, 0.05, put) // mark below intrinsic (call ITM / put OTM etc.)
			add(s, s*1.05, day, 0.05, put)
			add(s, s*0.5, 21.0/365.25, 0.05, put)
			add(s, s*1.5, 21.0/365.25, 400, put)
			add(s, s, 21.0/365.25, 400, put)
		}
	}
	// solver default start (estimate 0 -> 0.1) and fixed-vol delta edges
	zero := 0.0
	for _, put := range []bool{false, true} {
		for _, k := range []float64{7600, 7700, 7800} {
			addE(7700, k, 1.0/365.25, 5, put, 0, nil)
			addE(7700, k, 1.0/365.25, 50, put, 0, nil)
			addE(7700, k, 1.0/365.25, 5, put, 0.11, &zero) // zero vol: w is NaN or +-Inf
		}
	}
	out, _ := json.Marshal(rows)
	os.WriteFile(os.Args[1], out, 0644)
	greeks(os.Args[2])
}

type GRow struct {
	S, K, T, V, R       float64
	Put                 bool
	Scale               int
	GammaBits, VegaBits uint64
	ThetaBits, RhoBits  uint64
}

// greeks writes gamma/vega/theta/rho rows. v1.7.0 GetTheta/GetRho test "CALL" (v1.1.0 tested
// lowercase "call"); the Rust API takes a bool.
func greeks(path string) {
	rng := rand.New(rand.NewSource(43))
	var rows []GRow
	add := func(s, k, t, v float64, put bool, scale int) {
		const r = 0.0425
		ct := "CALL"
		if put {
			ct = "PUT"
		}
		rows = append(rows, GRow{s, k, t, v, r, put, scale,
			math.Float64bits(blackscholes.GetGamma(s, k, t, v, r)),
			math.Float64bits(blackscholes.GetVega(s, k, t, v, r)),
			math.Float64bits(blackscholes.GetTheta(s, k, t, v, r, ct, scale)),
			math.Float64bits(blackscholes.GetRho(s, k, t, v, r, ct, scale))})
	}
	scales := []int{0, 0, 252, 100, 1}
	for i := 0; i < 4000; i++ {
		s := 7700.0
		if i%2 == 1 {
			s = 770
		}
		s += (rng.Float64() - 0.5) * 0.02 * s
		k := math.Round(s*(1+(rng.Float64()*0.3-0.15))/5) * 5
		t := math.Exp(math.Log(1.0/(365.25*24*60)) + rng.Float64()*(math.Log(21.0/365.25)-math.Log(1.0/(365.25*24*60))))
		add(s, k, t, 0.02+rng.Float64()*0.8, i%4 >= 2, scales[i%5])
	}
	day := 1.0 / 365.25
	minute := day / 1440
	for _, put := range []bool{false, true} {
		for _, sc := range []int{0, 7} {
			for _, s := range []float64{7700, 770} {
				add(s, s, day, 0.2, put, sc)
				add(s, s, minute/60, 0.2, put, sc)
				add(s, s, 1e-12, 0.2, put, sc)
				add(s, s, day, 1e-9, put, sc)
				add(s, s*0.7, day, 0.2, put, sc)
				add(s, s*1.3, day, 0.2, put, sc)
				add(s, s*0.5, 21.0/365.25, 0.05, put, sc)
				add(s, s*1.5, 21.0/365.25, 0.05, put, sc)
				add(s, s, day, 0, put, sc)      // zero vol, k == s: w NaN
				add(s*1.01, s, day, 0, put, sc) // zero vol: w +Inf
				add(s*0.99, s, day, 0, put, sc) // zero vol: w -Inf
				add(s, s, 0, 0.2, put, sc)      // zero t, k == s: w NaN
				add(s, s*1.01, 0, 0.2, put, sc) // zero t: w Inf
			}
		}
	}
	out, _ := json.Marshal(rows)
	os.WriteFile(path, out, 0644)
}
