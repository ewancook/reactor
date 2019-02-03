package main

import (
	"math"
	"testing"
)

const tolerance = 0.01

func TestReaction1Enthalpy(t *testing.T) {
	res, expected := reaction1Enthalpy(298.15), 206.20
	if math.Abs(res-expected) >= tolerance {
		t.Errorf("incorrect reaction enthalpy: expected %f; got %f", expected, res)
	}
}

func TestReaction2Enthalpy(t *testing.T) {
	res, expected := reaction2Enthalpy(298.15), -41.20
	if math.Abs(res-expected) >= tolerance {
		t.Errorf("incorrect reaction enthalpy: expected %f; got %f", expected, res)
	}
}

func TestReaction3Enthalpy(t *testing.T) {
	res, expected := reaction3Enthalpy(298.15), 165.00
	if math.Abs(res-expected) >= tolerance {
		t.Errorf("incorrect reaction enthalpy: expected %f; got %f", expected, res)
	}
}
