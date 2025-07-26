import numpy as np

class MaterialBalance():

	@staticmethod
	def gas_solubility(rho,fvf,api,gg):
		"""
		The gas solubility can also be calculated rigorously from the experimental
		measured PVT data at the specified pressure and temperature.

		The method relates the gas solubility Rs to oil density, specific
		gravity of the oil, gas gravity, and the oil formation volume factor

		rho = oil density, lb/ft3
		fvf = oil formation volume factor, bbl/STB

		go = specific gravity of the stock-tank oil
		gg = specific gravity of the solution gas

		McCain (1991) pointed out that the weight average of separator and
		stock-tank gas specific gravities should be used for gg. The error in calculating
		Rs by using the above equation will depend only on the accuracy
		of the available PVT data.

		"""
		go = 141.5/(api+131.5)

		return (fvf*rho-62.4*go)/(0.0136*gg)

if __name__ == "__main__":

	print(GasSolubility.get(38.13,1.528,47.1,0.851))
	print(GasSolubility.get(40.95,1.474,40.7,0.855))
	print(GasSolubility.get(37.37,1.529,48.6,0.911))
	print(GasSolubility.get(38.92,1.619,40.5,0.898))
	print(GasSolubility.get(37.70,1.570,44.2,0.781))
	print(GasSolubility.get(46.79,1.385,27.3,0.848))