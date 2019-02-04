package main

import (
	"flag"
	"math"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
	"github.com/cpmech/gosl/plt"
)

func main() {
	flowCH4 := flag.Float64("CH4", 125, "initial flow of methane (mol/s)")
	flowH2 := flag.Float64("H2", 16.25, "initial flow of hydrogen (mol/s)")
	flowCO := flag.Float64("CO", 14.6, "initial flow of carbon monoxide (mol/s)")
	flowCO2 := flag.Float64("CO2", 2.71, "initial flow of carbon dioxide (mol/s)")
	flowH2O := flag.Float64("H2O", 140, "initial flow of steam (mol/s)")
	ρb := flag.Float64("density", 1.6, "bulk density (kg/m^3)")
	U := flag.Float64("U", 50, "heat transfer coefficient (W/Km^2)")
	T := flag.Float64("T", 1400, "initial reactor temperature (K)")
	Tα := flag.Float64("Talpha", 1500, "heating gas temperature, Tα (K)")
	P := flag.Float64("P", 2000, "initial reactor pressure (kPa)")
	D := flag.Float64("D", 5, "reactor diameter (m)")
	ϕ := flag.Float64("voidage", 0.4, "bed voidage (ϕ)")
	μ := flag.Float64("viscosity", 0.00008, "gas viscosity (μ)")
	Dp := flag.Float64("Dp", 0.010, "particle diameter (m)")
	ρc := flag.Float64("catalyst-density", 1140, "catalyst-density (kg/m^3)")
	W := flag.Float64("W", 100, "catalyst weight (kg)")
	flag.Parse()

	F0 := *flowCO + *flowH2 + *flowCH4 + *flowCO2 + *flowH2O

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
			partials[compound] = flow/totalFlow* y[6]
		}
		G := (16.04*flows["CH4"] + 44.01*flows["CO2"] + 28.01*flows["CO"] + 2.016*flows["H2"] + 18.015*flows["H2O"]) / 1000
		area := math.Pi*math.Pow(*D, 2)/4
		denominator := _denominator(y[5], partials)
		alpha := α(β(*ϕ, G, *Dp, *μ, *ρb), area, *ρc, *ϕ, *P)

		f[0] = dFCOdW(y[5], denominator, partials)
		f[1] = dFH2dW(y[5], denominator, partials)
		f[2] = dFCH4dW(y[5], denominator, partials)
		f[3] = dFCO2dW(y[5], denominator, partials)
		f[4] = dFH2OdW(y[5], denominator, partials)
		f[5] = dTdW(*U, *D, *ρb, *Tα, y[5], denominator, flows, partials)
		f[6] = dPdW(alpha, y[6], *P, y[5], *T, totalFlow, F0)
	}

	config := ode.NewConfig("dopri5", "", nil)
	config.SetStepOut(true, nil)

	parameters := []float64{*flowCO, *flowH2, *flowCH4, *flowCO2, *flowH2O, *T, *P}
	solver := ode.NewSolver(len(parameters), config, ODEs, nil, nil)
	defer solver.Free()
	solver.Solve(la.NewVectorSlice(parameters), 0, *W)

	wValues := solver.Out.GetStepX()
	yValues := solver.Out.GetStepYtableT()

	conversion := make([]float64, len(wValues))
	for i := 0; i < len(wValues); i++ {
		conversion[i] = 1 - yValues[2][i]/yValues[2][0]
	}
	plt.Subplot(2, 2, 1)
	plt.Plot(wValues,conversion, nil)
	plt.Grid(nil)
	plt.AxisRange(0, *W, 0, 1)
	plt.SetXlabel("Catalyst (kg)", nil)
	plt.SetYlabel("Methane Conversion", nil)

	plt.Subplot(2, 2, 2)
	plt.Plot(wValues, yValues[5], nil)
	plt.Grid(nil)
	plt.AxisXmax(*W)
	plt.SetXlabel("Catalyst (kg)", nil)
	plt.SetYlabel("Temperature (K)", nil)

	plt.Subplot(2, 2, 3)
	plt.Plot(wValues, yValues[6], nil)
	plt.Grid(nil)
	plt.AxisXmax(*W)
	plt.SetXlabel("Catalyst (kg)", nil)
	plt.SetYlabel("Pressure (kPa)", nil)
	plt.Show()
}
