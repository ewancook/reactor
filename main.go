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
	flowCH4 := flag.Float64("CH4", 106, "initial flow of methane (mol/s)")
	flowH2 := flag.Float64("H2", 6.57, "initial flow of hydrogen (mol/s)")
	flowCO := flag.Float64("CO", 0.001, "initial flow of carbon monoxide (mol/s)")
	flowCO2 := flag.Float64("CO2", 2.988, "initial flow of carbon dioxide (mol/s)")
	flowH2O := flag.Float64("H2O", 383, "initial flow of steam (mol/s)")
	flowC2H6 := flag.Float64("C2H6", 10, "initial flow of ethane (mol/s)")
	ρ := flag.Float64("density", 6.38, "gas density (kg/m^3)")
	U := flag.Float64("U", 40, "heat transfer coefficient (W/Km^2)")
	T := flag.Float64("T", 823.15, "initial reactor temperature (K)")
	Tα := flag.Float64("Talpha", 2000, "heating gas temperature, Tα (K)")
	P := flag.Float64("P", 2350, "initial reactor pressure (kPa)")
	D := flag.Float64("D", 0.11, "reactor diameter (m)")
	ϕ := flag.Float64("voidage", 0.44, "bed voidage (ϕ)")
	μ := flag.Float64("viscosity", 0.00002, "gas viscosity (μ)")
	Dp := flag.Float64("Dp", 0.013, "particle diameter (m)")
	ρb := flag.Float64("catalyst-density", 870, "catalyst-density (kg/m^3)")
	l := flag.Float64("l", 15, "tube length (m)")
	t := flag.Float64("tubes", 200, "number of tubes")
	nograph := flag.Bool("nograph", false, "stops plotting of graphs")

	// flue gases
	flueN2 := flag.Float64("flueN2", 738.5, "flue flowrate of nitrogen (mol/s)")
	flueCO2 := flag.Float64("flueCO2", 137.15, "flue flowrate of carbon dioxide (mol/s)")
	flueH2O := flag.Float64("flueH2", 137.15, "flue flowrate of steam (mol/s)")
	flueO2 := flag.Float64("flueCH4", 42.2, "flue flowrate of oxygen (mol/s)")

	flag.Parse()

	ρc := *ρb / (1.0 - *ϕ)

	F0 := *flowCO + *flowH2 + *flowCH4 + *flowCO2 + *flowH2O + *flowC2H6
	area := math.Pi * math.Pow(*D, 2) / 4
	W := *ρb * area * *l
	totalFlue := *flueN2 + *flueH2O + *flueCO2 + *flueO2
	ODEs := func(f la.Vector, h, x float64, y la.Vector) {
		totalFlow := y[0] + y[1] + y[2] + y[3] + y[4] + y[5]
		flows := map[string]float64{
			"CO":   y[0],
			"H2":   y[1],
			"CH4":  y[2],
			"CO2":  y[3],
			"H2O":  y[4],
			"C2H6": y[5],
		}
		partials := map[string]float64{}
		for compound, flow := range flows {
			partials[compound] = flow / totalFlow * y[7]
		}
		G := (16.04*flows["CH4"] + 44.01*flows["CO2"] + 28.01*flows["CO"] + 2.016*flows["H2"] + 18.015*flows["H2O"] + 30.069*flows["C2H6"]) / 1000 / area
		denominator := _denominator(y[6], partials)
		beta := β(*ϕ, G, *Dp, *μ, *ρ)
		alpha := α(beta, area, ρc, *ϕ, *P)

		aveCP := (thermo.SpecificHeat("N2", y[6])**flueN2 +
			thermo.SpecificHeat("CO2", y[6])**flueCO2 +
			thermo.SpecificHeat("H2O", y[6])**flueH2O +
			thermo.SpecificHeat("O2", y[6])**flueO2) / totalFlue

		f[0] = dFCOdW(y[6], denominator, partials)
		f[1] = dFH2dW(y[6], denominator, partials)
		f[2] = dFCH4dW(y[6], denominator, partials)
		f[3] = dFCO2dW(y[6], denominator, partials)
		f[4] = dFH2OdW(y[6], denominator, partials)
		f[5] = dFC2H6dW(y[6], partials)
		f[6] = dTdW(*U, *D, *ρb, y[8], y[6], denominator, flows, partials)
		f[7] = dPdW(alpha, y[7], *P, y[6], *T, totalFlow, F0)
		f[8] = dTαdW(*U, *D, *ρb, y[6], y[8], totalFlue, aveCP) * *t
	}

	config := ode.NewConfig("radau5", "", nil)
	config.SetStepOut(true, nil)

	parameters := []float64{*flowCO / *t,
		*flowH2 / *t,
		*flowCH4 / *t,
		*flowCO2 / *t,
		*flowH2O / *t,
		*flowC2H6 / *t, *T, *P, *Tα}
	F0 = F0 / *t
	solver := ode.NewSolver(len(parameters), config, ODEs, nil, nil)
	defer solver.Free()
	solver.Solve(la.NewVectorSlice(parameters), 0, W)

	wValues := solver.Out.GetStepX()
	yValues := solver.Out.GetStepYtableT()

	conversion := 1 - yValues[2][len(wValues)-1]/yValues[2][0]
	pressureDrop := *P - yValues[7][len(wValues)-1]

	fmt.Printf("tubes: %.0f; conversion: %.2f; pressure drop (kPa): %.4f; outlet temperature %2f (K)\n", *t, conversion, pressureDrop, yValues[6][len(wValues)-1])
	flows := []float64{parameters[0] * *t,
		parameters[1] * *t,
		parameters[2] * *t,
		parameters[3] * *t,
		parameters[4] * *t,
		parameters[5] * *t,
	}

	fmt.Printf("flows (mol/s); CO: %.2f; H2: %.2f; CH4: %.2f; CO2: %.2f; H2O %.2f; C2H6: %.2f\n", flows[0], flows[1], flows[2], flows[3], flows[4], flows[5])
	var conversions []float64
	var ethaneConversions []float64

	for _, v := range yValues[2] {
		conversions = append(conversions, 1-v/yValues[2][0])

	}

	for _, e := range yValues[5] {
		ethaneConversions = append(ethaneConversions, 1-e/yValues[5][0])
	}

	if *nograph {
		return
	}

	plt.Subplot(2, 3, 1)
	plt.Plot(wValues, conversions, nil)
	plt.Grid(nil)
	plt.SetLabels("Catalyst (kg)", "Methane Conversion", nil)

	plt.Subplot(2, 3, 2)
	plt.Plot(wValues, yValues[7], nil)
	plt.SetTicksNormal()
	plt.Grid(nil)
	plt.SetLabels("Catalyst (kg)", "Presssure (kPa)", nil)

	plt.Subplot(2, 3, 3)
	plt.Plot(wValues, yValues[8], nil)
	plt.Grid(nil)
	plt.SetLabels("Catalyst (kg)", "Talpha (K)", nil)

	plt.Subplot(2, 3, 4)
	plt.Plot(wValues, yValues[6], nil)
	plt.Grid(nil)
	plt.SetLabels("Catalyst (kg)", "T (K)", nil)

	plt.Subplot(2, 3, 5)
	plt.Plot(wValues, ethaneConversions, nil)
	plt.AxisYmin(0)
	plt.Grid(nil)
	plt.SetLabels("Catalyst (kg)", "Ethane Conversion", nil)

	plt.Show()
}
