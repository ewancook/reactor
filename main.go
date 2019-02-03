package main

import (
	"flag"
	"fmt"
	"math"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
	"github.com/ewancook/reactor/thermo"
)

func main() {
	flowCH4 := flag.Float64("CH4", 0, "initial flow of methane (mol/s)")
	flowH2 := flag.Float64("H2", 0, "initial flow of hydrogen (mol/s)")
	flowCO := flag.Float64("CO", 0, "initial flow of carbon monoxide (mol/s)")
	flowCO2 := flag.Float64("CO2", 0, "initial flow of carbon dioxide (mol/s)")
	flowH2O := flag.Float64("H2O", 0, "initial flow of steam (mol/s)")
	/*density := flag.Float64("density", 1000, "bulk density (kg/m^3)")
	heatTransfer := flag.Float64("U", 30, "heat transfer coefficient (W/Km^2)")
	temperature := flag.Float64("T", 298.15, "initial reactor temperature (K)")
	temperatureAlpha := flag.Float64("Talpha", 298.15, "heating gas temperature, Tα (K)")
	pressure := flag.Float64("P", 101.325, "initial reactor pressure (kPa)")
	*/
	flag.Parse()

	totalFlow := *flowCO + *flowH2 + *flowCH4 + *flowCO2 + *flowH2O
	fmt.Println(totalFlow)
	y := la.NewVectorSlice([]float64{4})
	sol := ode.NewSolver(len(y), ode.NewConfig("dopri5", "", nil), solveThis, nil, nil)
	defer sol.Free()
	sol.Solve(y, 0, 1)
	fmt.Println(β(0.4, 50.0/3600, 0.005, 2.3*math.Pow(10, -5), 12))

	fmt.Println(reaction1Enthalpy(298.15))
	fmt.Println(reaction2Enthalpy(298.15))
	fmt.Println(reaction3Enthalpy(298.15))
}

func dTdW(U, D, ρb, Tα, T float64, flows map[string]float64) float64 {
	var denominator float64
	for compound, flow := range flows {
		denominator += thermo.SpecificHeat(compound, T) * flow
	}
	return 0
}

func dPdW(α, P, initialP, T, initialT, F, initialF float64) float64 {
	return -α / 2 * initialP / (P / initialP) * (T / initialT) * (F / initialF)
}

func α(β, area, ρ, ϕ, P float64) float64 {
	return 2 * β / (area * ρ * (1 - ϕ) * P)
}

func β(ϕ, G, Dp, μ, ρ float64) float64 {
	return (G * (1 - ϕ) / (ρ * Dp * math.Pow(ϕ, 3))) * (1.75*G + 150*(1-ϕ)*μ/Dp)
}

func solveThis(f la.Vector, h, x float64, y la.Vector) {
	f[0] = math.Sin(x+y[0]) - math.Exp(x)
}
