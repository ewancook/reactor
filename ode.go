package main

import (
	"math"

	"github.com/ewancook/reactor/thermo"
)

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
	return (U*(4/D)/ρb*(Tα-T) - heats*1000) / denominator
}

func dPdW(alpha, P, P0, T, T0, F, F0 float64) float64 {
	return -alpha / 2 * P0 / (P / P0) * (T / T0) * (F / F0)
}

func α(beta, area, ρc, ϕ, P0 float64) float64 {
	return 2 * beta / (area * ρc * (1 - ϕ) * P0 * 1000)
}

func β(ϕ, G, Dp, μ, ρg float64) float64 {
	return (G * (1 - ϕ) / (ρg * Dp * math.Pow(ϕ, 3))) * (1.75*G + 150*(1-ϕ)*μ/Dp)
}

func dTαdW(U, D, ρb, T, Tα, mc, aveCP float64) float64 {
	return U * 4 / D / ρb * (T - Tα) / (mc * aveCP)
}
