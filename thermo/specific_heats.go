package thermo

import (
	"math"
)

var cpData = map[string]func(float64) float64{
	"CO":   cpCarbonMonoxide,
	"H2O":  cpSteam,
	"H2":   cpHydrogen,
	"CO2":  cpCarbonDioxide,
	"CH4":  cpMethane,
	"N2":   cpNitrogen,
	"O2":   cpOxygen,
	"C2H6": cpEthane,
}

func SpecificHeat(compound string, T float64) float64 {
	return cpData[compound](T)
}

func cpIntegral(T, A, B, C, D, E float64) float64 {
	T = T / 1000
	return A + B*T + C*math.Pow(T, 2) + D*math.Pow(T, 3) + E/math.Pow(T, 2)
}

func cpCarbonMonoxide(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateCarbonMonoxide(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpSteam(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateSteam(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpHydrogen(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateHydrogen(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpCarbonDioxide(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateCarbonDioxide(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpMethane(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateMethane(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpNitrogen(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateNitrogen(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpOxygen(T float64) float64 {
	a, b, c, d, e, _, _, _ := shomateOxygen(T)
	return cpIntegral(T, a, b, c, d, e)
}

func cpEthane(T float64) float64 {
	return 4*math.Pow(10, -18)*math.Pow(T, 6) -
		4*math.Pow(10, -14)*math.Pow(T, 5) + 2*math.Pow(10, -10)*math.Pow(T, 4) -
		3*math.Pow(10, -7)*math.Pow(T, 3) + 0.0003*math.Pow(T, 2) + 0.0308*T + 29.265
}
