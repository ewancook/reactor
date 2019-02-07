package thermo

import (
	"math"
	"testing"
)

const tolerance = 0.05

func _compareSpecificHeat(t *testing.T, res, expected float64) {
	if math.Abs(res-expected) >= tolerance {
		t.Errorf("incorrect specific heat: expected %f; got %f", expected, res)
	}
}

func TestCarbonMonoxideSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		298:  29.15,
		600:  30.47,
		1000: 33.18,
		2000: 36.21,
		6000: 38.37,
	}
	for T, expected := range results {
		res := SpecificHeat("CO", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestSteamSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		800:  38.74,
		1300: 44.94,
		1700: 48.92,
		2000: 51.20,
		3000: 55.74,
		6000: 60.59,
	}
	for T, expected := range results {
		res := SpecificHeat("H2O", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestHydrogenSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		298:  28.84,
		600:  29.32,
		1000: 30.20,
		1500: 32.30,
		2000: 34.28,
		6000: 41.97,
	}
	for T, expected := range results {
		res := SpecificHeat("H2", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestCarbonDioxideSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		298:  37.12,
		600:  47.32,
		1000: 54.30,
		2000: 60.34,
		4000: 63.25,
		6000: 64.98,
	}
	for T, expected := range results {
		res := SpecificHeat("CO2", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestMethaneSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		298:  35.64,
		600:  52.23,
		1000: 71.79,
		2000: 94.38,
		4000: 104.2,
		6000: 106.4,
	}
	for T, expected := range results {
		res := SpecificHeat("CH4", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestNitrogenSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		200:  29.11,
		500:  29.58,
		1000: 32.69,
		2000: 35.97,
		4000: 37.55,
		6000: 38.27,
	}
	for T, expected := range results {
		res := SpecificHeat("N2", T)
		_compareSpecificHeat(t, res, expected)
	}
}

func TestOxygenSpecificHeat(t *testing.T) {
	results := map[float64]float64{
		200:  29.11,
		600:  32.11,
		1000: 34.86,
		2000: 37.75,
		4000: 41.42,
		6000: 44.39,
	}
	for T, expected := range results {
		res := SpecificHeat("O2", T)
		_compareSpecificHeat(t, res, expected)
	}
}
