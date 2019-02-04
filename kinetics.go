package main

import "math"

const R = 8.314

func k1(T float64) float64 {
	return 5.922 * math.Pow(10, 8) * math.Exp(-209200/R/T)
}

func k2(T float64) float64 {
	return 6.028 * math.Pow(10, -4) * math.Exp(-15400/R/T)
}

func k3(T float64) float64 {
	return 1.093 * math.Pow(10, 3) * math.Exp(-109400/R/T)
}

func kCO(T float64) float64 {
	return 5.127 * math.Pow(10, -13) * math.Exp(140000/R/T)
}

func kH2(T float64) float64 {
	return 5.68 * math.Pow(10, -10) * math.Exp(93400/R/T)
}

func kH2O(T float64) float64 {
	return 9.251 * math.Exp(-15900/R/T)
}

func kp1(T float64) float64 {
	return 1.2 * math.Pow(10, 17) * math.Exp(-26830/T)
}

func kp2(T float64) float64 {
	return 1.8 * math.Pow(10, -2) * math.Exp(4400/T)
}

func kp3(T float64) float64 {
	return 2.1 * math.Pow(10, 15) * math.Exp(-22430/T)
}
