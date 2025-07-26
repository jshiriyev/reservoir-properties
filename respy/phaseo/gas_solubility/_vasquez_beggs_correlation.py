import numpy as np

def vasquez_beggs(T,p,API,gg,Tsep,psep):
	"""
	Vasquez and Beggs (1980) presented an improved empirical correlation
	for estimating Rs. The correlation was obtained by regression analysis
	using 5,008 measured gas solubility data points.

	Based on oil gravity, the measured data were divided into two groups.
	This division was made at a value of oil gravity of 30°API.

	Realizing that the value of the specific gravity of the gas depends on
	the conditions under which it is separated from the oil, Vasquez and
	Beggs proposed that the value of the gas specific gravity as obtained
	from a separator pressure of 100 psig be used in the above equation. This
	reference pressure was chosen because it represents the average field
	separator conditions.

	ggs = gas gravity at the reference separator pressure
	gg = gas gravity at the actual separator conditions of psep and Tsep
	psep = actual separator pressure, psia
	Tsep = actual separator temperature, °F

	"""
	if API<30:
		C1 = 0.0362
		C2 = 1.0937
		C3 = 25.7240
	else:
		C1 = 0.0178
		C2 = 1.1870
		C3 = 23.931

	ggs = gg*(1+5.912e-5*API*Tsep*np.log10((psep+14.7)/114.7))

	print(f"\n{ggs=}")

	return C1*ggs*(p+14.7)**C2*np.exp(C3*API/(T+460))

if __name__ == "__main__":

	print(vasquez_beggs(250,2377,47.1,0.851, 60,150))
	print(vasquez_beggs(220,2620,40.7,0.855, 75,100))
	print(vasquez_beggs(260,2051,48.6,0.911, 72,100))
	print(vasquez_beggs(237,2884,40.5,0.898,120, 60))
	print(vasquez_beggs(218,3045,44.2,0.781, 60,200))
	print(vasquez_beggs(180,4239,27.3,0.848,173, 85))