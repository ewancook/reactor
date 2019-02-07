package thermo

func shomateCarbonMonoxide(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 1300:
		return 25.56759, 6.09613, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665, -110.5271
	default:
		return 35.15070, 1.300095, -0.205921, 0.013550, -3.282780, -127.8375, 231.7120, -110.5271
	}
}

func shomateSteam(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 1700:
		return 30.092, 6.832514, 6.793435, -2.53448, 0.082139, -250.881, 223.3967, -241.8264
	default:
		return 41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264
	}
}

func shomateHydrogen(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 1000:
		return 33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0
	case 1000 < T && T <= 2500:
		return 18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0
	default:
		return 43.413560, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354, 0.0
	}
}

func shomateCarbonDioxide(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 1200:
		return 24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224
	default:
		return 58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224
	}
}

func shomateMethane(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 1300:
		return -0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376, 158.7163, -74.87310
	default:
		return 85.81217, 11.26467, -2.114146, 0.138190, -26.42221, -153.5327, 224.4143, -74.87310
	}
}

func shomateNitrogen(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 500:
		return 28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0
	case T > 500 && T <= 2000:
		return 19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0
	default:
		return 35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810, 0.0
	}
}

func shomateOxygen(T float64) (float64, float64, float64, float64, float64, float64, float64, float64) {
	switch {
	case T <= 700:
		return 31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0
	case T > 700 && T <= 2000:
		return 30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0
	default:
		return 20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0
	}
}
