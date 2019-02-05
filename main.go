package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
	"github.com/cpmech/gosl/plt"
	"github.com/ewancook/reactor/thermo"
)

func main() {
	flowCH4 := flag.Float64("CH4", 125, "initial flow of methane (mol/s)")
	flowH2 := flag.Float64("H2", 16.25, "initial flow of hydrogen (mol/s)")
	flowCO := flag.Float64("CO", 14.6, "initial flow of carbon monoxide (mol/s)")
	flowCO2 := flag.Float64("CO2", 2.71, "initial flow of carbon dioxide (mol/s)")
	flowH2O := flag.Float64("H2O", 383, "initial flow of steam (mol/s)")
	ρ := flag.Float64("density", 4, "gas density (kg/m^3)")
	U := flag.Float64("U", 50, "heat transfer coefficient (W/Km^2)")
	T := flag.Float64("T", 630, "initial reactor temperature (K)")
	Tα := flag.Float64("Talpha", 2000, "heating gas temperature, Tα (K)")
	P := flag.Float64("P", 2500, "initial reactor pressure (kPa)")
	D := flag.Float64("D", 0.11, "reactor diameter (m)")
	ϕ := flag.Float64("voidage", 0.44, "bed voidage (ϕ)")
	μ := flag.Float64("viscosity", 0.00002, "gas viscosity (μ)")
	Dp := flag.Float64("Dp", 0.010, "particle diameter (m)")
	ρc := flag.Float64("catalyst-density", 1140, "catalyst-density (kg/m^3)")
	l := flag.Float64("l", 15, "tube length (m)")
	flag.Parse()

	ρb := (1.0 - *ϕ) * *ρc

	F0 := *flowCO + *flowH2 + *flowCH4 + *flowCO2 + *flowH2O
	area := math.Pi * math.Pow(*D, 2) / 4
	W := ρb * area * *l

	ODEs := func(f la.Vector, h, x float64, y la.Vector) {
		totalFlow := y[0] + y[1] + y[2] + y[3] + y[4]
		flows := map[string]float64{
			"CO":  y[0],
			"H2":  y[1],
			"CH4": y[2],
			"CO2": y[3],
			"H2O": y[4],
		}
		partials := map[string]float64{}
		for compound, flow := range flows {
			partials[compound] = flow / totalFlow * y[6]
		}
		G := (16.04*flows["CH4"] + 44.01*flows["CO2"] + 28.01*flows["CO"] + 2.016*flows["H2"] + 18.015*flows["H2O"]) / 1000
		denominator := _denominator(y[5], partials)
		alpha := α(β(*ϕ, G, *Dp, *μ, *ρ), area, *ρc, *ϕ, *P)

		var totalCP float64
		for _, c := range []string{"CO", "CO2", "H2", "H2O", "CH4"} {
			totalCP += thermo.SpecificHeat(c, y[5])
		}
		aveCP := totalCP / 5

		f[0] = dFCOdW(y[5], denominator, partials)
		f[1] = dFH2dW(y[5], denominator, partials)
		f[2] = dFCH4dW(y[5], denominator, partials)
		f[3] = dFCO2dW(y[5], denominator, partials)
		f[4] = dFH2OdW(y[5], denominator, partials)
		f[5] = dTdW(*U, *D, ρb, y[7], y[5], denominator, flows, partials)
		f[6] = dPdW(alpha, y[6], *P, y[5], *T, totalFlow, F0)
		f[7] = dTαdW(*U, *D, ρb, y[5], y[7], totalFlow, aveCP)
	}

	config := ode.NewConfig("dopri5", "", nil)
	config.SetStepOut(true, nil)

	// iteration for 'best'
	tubes := []float64{}
	conversions := []float64{}
	pressures := []float64{}

	for t := 50.0; t <= 600.0; t++ {
		parameters := []float64{*flowCO / t,
			*flowH2 / t,
			*flowCH4 / t,
			*flowCO2 / t,
			*flowH2O / t, *T, *P, *Tα}
		solver := ode.NewSolver(len(parameters), config, ODEs, nil, nil)
		defer solver.Free()
		solver.Solve(la.NewVectorSlice(parameters), 0, W)

		wValues := solver.Out.GetStepX()
		yValues := solver.Out.GetStepYtableT()

		conversion := 1 - yValues[2][len(wValues)-1]/yValues[2][0]
		pressureDrop := yValues[6][0] - yValues[6][len(wValues)-1]
		tubes = append(tubes, t)
		conversions = append(conversions, conversion)
		pressures = append(pressures, pressureDrop)
		if conversion >= 0.75 {
			fmt.Printf("tubes: %.0f; pressure drop (kPa): %.2e\n", t, pressureDrop)
			break
		}
	}
	plt.Subplot(1, 2, 1)
	plt.Plot(tubes, conversions, nil)
	plt.Grid(nil)
	plt.SetLabels("Number of Tubes", "Conversion", nil)

	plt.Subplot(1, 2, 2)
	plt.Plot(tubes, pressures, nil)
	plt.Grid(nil)
	plt.SetLabels("Number of Tubes", "Presssure Drop (kPa)", nil)
	plt.Show()
}
