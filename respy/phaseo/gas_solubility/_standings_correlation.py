def standing(T,p,API,gg):
	"""
	Standing (1947) proposed a graphical correlation for determining the
	gas solubility as a function of pressure, gas specific gravity, API gravity,
	and system temperature.

	The correlation was developed from a total of 105 experimentally determined
	data points on 22 hydrocarbon mixtures from California crude oils and
	natural gases. The proposed correlation has an average error of 4.8%.
	
	Standing (1981) expressed his proposed graphical correlation in the following
	more convenient mathematical form:

	T = temperature, °F
	p = system pressure, psia

	gg = solution gas specific gravity

	It should be noted that Standing’s equation is valid for applications at
	and below the bubble-point pressure of the crude oil.

	"""
	x = 0.0125 * API - 0.00091 * T
	# print(f"\n{x=}")

	Rs = gg*((p/18.2+1.4)*10**x)**1.20482

	return Rs

if __name__ == "__main__":

	print(standing(250,2377+14.7,47.1,0.851))
	print(standing(220,2620+14.7,40.7,0.855))
	print(standing(260,2051+14.7,48.6,0.911))
	print(standing(237,2884+14.7,40.5,0.898))
	print(standing(218,3045+14.7,44.2,0.781))
	print(standing(180,4239+14.7,27.3,0.848))