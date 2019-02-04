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
	flowCH4 := flag.Float64("CH4", 125, "initial flow of methane (mol/s)")
	flowH2 := flag.Float64("H2", 16.25, "initial flow of hydrogen (mol/s)")
	flowCO := flag.Float64("CO", 14.6, "initial flow of carbon monoxide (mol/s)")
	flowCO2 := flag.Float64("CO2", 2.71, "initial flow of carbon dioxide (mol/s)")
	flowH2O := flag.Float64("H2O", 140, "initial flow of steam (mol/s)")
	ρb := flag.Float64("density", 1.6, "bulk density (kg/m^3)")
	U := flag.Float64("U", 30, "heat transfer coefficient (W/Km^2)")
	T := flag.Float64("T", 1133, "initial reactor temperature (K)")
	Tα := flag.Float64("Talpha", 1700, "heating gas temperature, Tα (K)")
	P := flag.Float64("P", 2000, "initial reactor pressure (kPa)")
	D := flag.Float64("D", 2, "reactor diameter (m)")
	ϕ := flag.Float64("voidage", 0.4, "bed voidage (ϕ)")
	μ := flag.Float64("viscosity", 0.00008, "gas viscosity (μ)")
	Dp := flag.Float64("Dp", 0.010, "particle diameter (m)")
	area := flag.Float64("area", .447, "cross-sectional area (m^2)")
	ρc := flag.Float64("catalyst-density", 1140, "catalyst-density (kg/m^3)")
	W := flag.Float64("W", 500, "catalyst weight (kg)")
	flag.Parse()
	G := (16.04**flowCH4 + 44.01**flowCO2 + 28.01**flowCO + 2.016**flowH2 + 18.015**flowH2O) / 1000

	P0 := *P
	T0 := *T
	F0 := *flowCO + *flowH2 + *flowCH4 + *flowCO2 + *flowH2O
	parameters := []float64{*flowCO, *flowH2, *flowCH4, *flowCO2, *flowH2O, *T, *P}
	fmt.Printf("%+v\n", parameters)

	solver := func(f la.Vector, h, x float64, y la.Vector) {
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
		denominator := _denominator(y[5], partials)
		alpha := α(β(*ϕ, G, *Dp, *μ, *ρb), *area, *ρc, *ϕ, P0)
		f[0] = dFCOdW(y[5], denominator, partials)
		f[1] = dFH2dW(y[5], denominator, partials)
		f[2] = dFCH4dW(y[5], denominator, partials)
		f[3] = dFCO2dW(y[5], denominator, partials)
		f[4] = dFH2OdW(y[5], denominator, partials)
		f[5] = dTdW(*U, *D, *ρb, *Tα, y[5], denominator, flows, partials)
		f[6] = dPdW(alpha, y[6], P0, y[5], T0, totalFlow, F0)
		var _, _, _, _, _, _ = T0, F0, alpha, U, Tα, D
		fmt.Printf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f \n", y[0], y[1], y[2], y[3], y[4], y[5], y[6]) //, P0, T0, F0, alpha, U, Tα, D)
	}
	config := ode.NewConfig("dopri5", "", nil)
	config.NmaxIt = 1e3
	sol := ode.NewSolver(len(parameters), config, solver, nil, nil)
	defer sol.Free()
	sol.Solve(la.NewVectorSlice(parameters), 0, *W)
	fmt.Printf("%+v\n", parameters)
}

func dFCH4dW(T, denominator float64, partials map[string]float64) float64 {
	return -reaction1(T, denominator, partials) - reaction3(T, denominator, partials)
}

func dFH2OdW(T, denominator float64, partials map[string]float64) float64 {
	return -reaction1(T, denominator, partials) -
		reaction2(T, denominator, partials) -
		2*reaction3(T, denominator, partials)
}

func dFH2dW(T, denominator float64, partials map[string]float64) float64 {
	return 3*reaction1(T, denominator, partials) +
		reaction2(T, denominator, partials) +
		4*reaction3(T, denominator, partials)
}

func dFCOdW(T, denominator float64, partials map[string]float64) float64 {
	return reaction1(T, denominator, partials) -
		reaction2(T, denominator, partials)
}

func dFCO2dW(T, denominator float64, partials map[string]float64) float64 {
	return reaction2(T, denominator, partials) +
		reaction3(T, denominator, partials)
}

func dTdW(U, D, ρb, Tα, T, reaction_denominator float64, flows, partials map[string]float64) float64 {
	var denominator float64
	for compound, flow := range flows {
		denominator += thermo.SpecificHeat(compound, T) * flow
	}
	heats := reaction1(T, reaction_denominator, partials)*reaction1Enthalpy(T) +
		reaction2(T, reaction_denominator, partials)*reaction2Enthalpy(T) +
		reaction3(T, reaction_denominator, partials)*reaction3Enthalpy(T)
	return (U*(4/D)/ρb*(Tα-T) - heats) / denominator
}

func dPdW(alpha, P, P0, T, T0, F, F0 float64) float64 {
	return -alpha / 2 * P0 / (P / P0) * (T / T0) * (F / F0)
}

func α(beta, area, ρc, ϕ, P0 float64) float64 {
	return 2 * beta / (area * ρc * (1 - ϕ) * P0)/102
}

func β(ϕ, G, Dp, μ, ρg float64) float64 {
	return (G * (1 - ϕ) / (ρg * Dp * math.Pow(ϕ, 3))) * (1.75*G + 150*(1-ϕ)*μ/Dp)
}
