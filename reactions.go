package main

import (
	. "math"

	. "github.com/ewancook/reactor/thermo"
)

func reaction1Enthalpy(T float64) float64 {
	return Enthalpy("CO", T) + Enthalpy("H2", T)*3 - Enthalpy("H2O", T) - Enthalpy("CH4", T)
}

func reaction2Enthalpy(T float64) float64 {
	return Enthalpy("H2", T) + Enthalpy("CO2", T) - Enthalpy("CO", T) - Enthalpy("H2O", T)
}

func reaction3Enthalpy(T float64) float64 {
	return Enthalpy("CO2", T) + Enthalpy("H2", T)*4 - Enthalpy("H2O", T)*2 - Enthalpy("CH4", T)
}

func reaction4Enthalpy(T float64) float64 {
	return Enthalpy("CO", T)*2 + Enthalpy("H2", T)*5 - Enthalpy("C2H6", T) - Enthalpy("H2O", T)*2
}

func _denominator(T float64, partials map[string]float64) float64 {
	return Pow(1+kCO(T)*partials["CO"]+kH2(T)*Pow(partials["H2"], 0.5)+kH2O(T)*partials["H2O"]/partials["H2"], 2)
}

func reaction1(T, denominator float64, partials map[string]float64) float64 {
	return (k1(T) * partials["CH4"] * Pow(partials["H2O"], 0.5) / Pow(partials["H2"], 1.25)) * (1 - partials["CO"]*Pow(partials["H2"], 3)/kp1(T)/partials["CH4"]/partials["H2O"]) / denominator
}

func reaction2(T, denominator float64, partials map[string]float64) float64 {
	return (k2(T) * partials["CO"] * Pow(partials["H2O"], 0.5) / Pow(partials["H2"], 0.5)) * (1 - partials["CO2"]*partials["H2"]/kp2(T)/partials["CO"]/partials["H2O"]) / denominator
}

func reaction3(T, denominator float64, partials map[string]float64) float64 {
	return (k3(T) * partials["CH4"] * partials["H2O"] / Pow(partials["H2"], 1.75)) * (1 - partials["CO2"]*Pow(partials["H2"], 4)/kp3(T)/partials["CH4"]/Pow(partials["H2O"], 2)) / denominator
}

func reaction4(T float64, partials map[string]float64) float64 {
	return (k4(T) * partials["C2H6"]) / Pow(1+25.2*partials["C2H6"]*partials["H2"]/partials["H2O"]+0.077*partials["H2O"]/partials["H2"], 2) / 3.6
}
