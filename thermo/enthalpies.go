package thermo

import (
	"math"
)

var enthalpyData = map[string]func(float64) float64{
	"CO":  hCarbonMonoxide,
	"H2O": hSteam,
	"H2":  hHydrogen,
	"CO2": hCarbonDioxide,
	"CH4": hMethane,
}

func Enthalpy(compound string, T float64) float64 {
	return enthalpyData[compound](T)
}

func hIntegral(T, A, B, C, D, E, F, G, H float64) float64 {
	T = T / 1000
	return A*T + math.Pow(T, 2)*B/2 + math.Pow(T, 3)*C/3 + math.Pow(T, 4)*D/4 - E/T + F - H
}

func hCarbonMonoxide(T float64) float64 {
	a, b, c, d, e, f, g, h := shomateCarbonMonoxide(T)
	return -110.5 + hIntegral(T, a, b, c, d, e, f, g, h)
}

func hSteam(T float64) float64 {
	a, b, c, d, e, f, g, h := shomateSteam(T)
	return -241.83 + hIntegral(T, a, b, c, d, e, f, g, h)
}

func hHydrogen(T float64) float64 {
	a, b, c, d, e, f, g, h := shomateHydrogen(T)
	return 0 + hIntegral(T, a, b, c, d, e, f, g, h)
}

func hCarbonDioxide(T float64) float64 {
	a, b, c, d, e, f, g, h := shomateCarbonDioxide(T)
	return -393.52 + hIntegral(T, a, b, c, d, e, f, g, h)
}

func hMethane(T float64) float64 {
	a, b, c, d, e, f, g, h := shomateMethane(T)
	return -74.87 + hIntegral(T, a, b, c, d, e, f, g, h)
}
