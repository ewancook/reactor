package thermo

import (
	"math"
	"testing"
)

func _compareEnthalpy(t *testing.T, res, expected float64) {
	if math.Abs(res-expected) >= tolerance {
		t.Errorf("incorrect enthalpy: expected %f; got %f", expected, res)
	}
}

func TestCarbonMonoxideEnthalpy(t *testing.T) {
	results := map[float64]float64{
		298:  0,
		600:  8.94,
		1000: 21.69,
		2000: 56.74,
		6000: 207.1,
	}
	for T, expected := range results {
		res := Enthalpy("CO", T)
		_compareEnthalpy(t, res+110.5, expected)
	}
}

func TestSteamEnthalpy(t *testing.T) {
	results := map[float64]float64{
		800:  18.00,
		1300: 38.94,
		1700: 57.76,
		2000: 72.79,
		3000: 126.5,
		6000: 302.3,
	}
	for T, expected := range results {
		res := Enthalpy("H2O", T)
		_compareEnthalpy(t, res+241.83, expected)
	}
}

func TestHydrogenEnthalpy(t *testing.T) {
	results := map[float64]float64{
		298:  0,
		600:  8.81,
		1000: 20.68,
		1500: 36.29,
		2000: 52.95,
		6000: 208.3,
	}
	for T, expected := range results {
		res := Enthalpy("H2", T)
		_compareEnthalpy(t, res, expected)
	}
}

func TestCarbonDioxideEnthalpy(t *testing.T) {
	results := map[float64]float64{
		298:  -0.01,
		600:  12.91,
		1000: 33.40,
		2000: 91.44,
		4000: 215.6,
		6000: 343.8,
	}
	for T, expected := range results {
		res := Enthalpy("CO2", T)
		_compareEnthalpy(t, res+393.52, expected)
	}
}

func TestMethaneEnthalpy(t *testing.T) {
	results := map[float64]float64{
		298:  0.00,
		600:  13.13,
		1000: 38.18,
		2000: 123.6,
		4000: 325.1,
		6000: 535.9,
	}
	for T, expected := range results {
		res := Enthalpy("CH4", T)
		_compareEnthalpy(t, res+74.87, expected)
	}
}
