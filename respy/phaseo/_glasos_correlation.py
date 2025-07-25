import numpy as np

class Glaso():

	def gas_solubility(T,p,api,gg):
		"""
		Glaso (1980) proposed a correlation for estimating the gas solubility as
		a function of the API gravity, pressure, temperature, and gas specific gravity.
		The correlation was developed from studying 45 North Sea crude oil
		samples. Glaso reported an average error of 1.28% with a standard deviation
		of 6.98%.

		T = temperature, °F
		p = system pressure, psia

		gg = solution gas specific gravity

		"""
		x = 2.8869 - (14.1811 - 3.3093*np.log10(p))**0.5
		# print(f"\n{x=}")

		pb_star = 10**x

		return gg*(api**0.989/T**0.172*pb_star)**1.2255

if __name__ == "__main__":

	print(glaso(250,2377+14.7,47.1,0.851))
	print(glaso(220,2620+14.7,40.7,0.855))
	print(glaso(260,2051+14.7,48.6,0.911))
	print(glaso(237,2884+14.7,40.5,0.898))
	print(glaso(218,3045+14.7,44.2,0.781))
	print(glaso(180,4239+14.7,27.3,0.848))