package main

import (
	. "github.com/ewancook/reactor/thermo"
)

func reaction1Enthalpy(T float64) (totalEnthalpy float64) {
	return Enthalpy("CO", T) + Enthalpy("H2", T)*3 - Enthalpy("H2O", T) - Enthalpy("CH4", T)
}

func reaction2Enthalpy(T float64) (totalEnthalpy float64) {
	return Enthalpy("H2", T) + Enthalpy("CO2", T) - Enthalpy("CO", T) - Enthalpy("H2O", T)
}

func reaction3Enthalpy(T float64) (totalEnthalpy float64) {
	return Enthalpy("CO2", T) + Enthalpy("H2", T)*4 - Enthalpy("H2O", T)*2 - Enthalpy("CH4", T)
}
