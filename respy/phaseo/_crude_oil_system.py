import numpy as np

class CrudeOilSystem:
	"""
	Petroleum (an equivalent term is crude oil) is a complex mixture
	consisting predominantly of hydrocarbons and containing sulfur, nitrogen,
	oxygen, and helium as minor constituents.

	The physical and chemical properties of crude oils vary considerably and
	are dependent on the concentration of the various types of hydrocarbons
	and minor constituents present.

	"""
	@staticmethod
	def get_sgco(rhoo:float|np.ndarray,rhow:float=62.4):
		"""
		The specific gravity of a crude oil (sgco) is defined as the ratio
		of the density of the oil to that of water. Both densities are measured
		at 60°F and atmospheric pressure:

		Inputs:
		------
		rhoo = Density of the crude oil, lb/ft3
		rhow = Density of the water, lb/ft3

		Returns:
		-------
		Specific gravity of the oil

		It should be pointed out that the liquid specific gravity is dimensionless,
		but traditionally is given the units 60°/60° to emphasize the fact that
		both densities are measured at standard conditions. The density of the
		water is approximately 62.4 lb/ft3.

		"""
		return rhoo/rhow

	@staticmethod
	def sgco_to_gAPI(sgco:float|np.ndarray):
		"""
		Calculates API gravity from specific gravity:
		
		Inputs:
		------
		sgco	: specific gravity of the oil

		Returns:
		-------
		gAPI 	: API gravity of the oil

		The API gravities of crude oils usually range from 47° API for the
		lighter crude oils to 10° API for the heavier asphaltic crude oils.

		"""
		return 141.5/sgco-131.5

	@staticmethod
	def gAPI_to_sgco(gAPI:float|np.ndarray):
		"""
		Calculates specific gravity from API gravity:
		
		Inputs:
		------
		gAPI 	: API gravity of the oil

		Returns:
		-------
		sgco 	: specific gravity of the oil

		"""
		return 141.5/(gAPI+131.5)

	@staticmethod
	def get_sgsg(*separators,ST:tuple[float,float]):
		"""
		Calculates the specific gravity of the solution gas.

		The specific gravity of the solution gas (sgsg) is described by the weighted
		average of the specific gravities of the separated gas from each separator.

		This weighted-average approach is based on the separator gas-oil ratio.
		
		The ST (stock-tank) is a tuple with Rst and gst values:
			Rst = gas-oil ratio from the stock tank, scf/ STB
			gst = gas gravity from the stock tank
		
		For each separator Rsep and gsep values must be provided:
			Rsep = separator gas-oil ratio, scf/STB
			gsep = separator gas gravity
		
		Returns:
		-------
		Specific gravity of the solution gas.

		"""
		Rst,gst = ST

		upper = Rst*gst
		lower = Rst

		for sep in separators:
			Rsep,gsep = sep

			upper += Rsep*gsep
			lower += Rsep

		return upper/lower

	def get_gassb():
		"""Calculates the stock-tank GOR. It should not be used
		if the separator temperature is >140F."""

		# CREEGER

		A1 =  0.3818
		A2 = -5.506
		A3 =  2.902
		A4 =  1.327
		A5 = -0.7355

		logrst = A1+A2*np.log(gamma_oil)+A3*np.log(gamma_gas_sp)+\
		            A4*np.log(Psp)+A5*np.log(Tsp)

		return np.exp(logrst)

	@staticmethod
	def get_gass(rho,fvf,sgsg,gAPI):
		"""
		The gas solubility can also be calculated rigorously from the experimental
		measured PVT data at the specified pressure and temperature.

		The method relates the gas solubility Rs to oil density, specific
		gravity of the oil, gas gravity, and the oil formation volume factor

		rho = oil density, lb/ft3
		fvf = oil formation volume factor, bbl/STB

		sgsg = specific gravity of the solution gas
		sgco = specific gravity of the stock-tank oil

		McCain (1991) pointed out that the weight average of separator and
		stock-tank gas specific gravities should be used for sgsg. The error in calculating
		Rs by using the above equation will depend only on the accuracy
		of the available PVT data.

		"""
		sgco = 141.5/(gAPI+131.5)

		return (fvf*rho-62.4*sgco)/(0.0136*sgsg)

if __name__ == "__main__":

	sgco = CrudeOilSystem.get_sgco(53)
	print(sgco)

	print(CrudeOilSystem.sgco_to_gAPI(sgco))

	print(CrudeOilSystem.get_sgsg((724,0.743),(202,0.956),ST=(58,1.296)))

	print(CrudeOilSystem.gass(38.13,1.528,47.1,0.851))
	print(CrudeOilSystem.gass(40.95,1.474,40.7,0.855))
	print(CrudeOilSystem.gass(37.37,1.529,48.6,0.911))
	print(CrudeOilSystem.gass(38.92,1.619,40.5,0.898))
	print(CrudeOilSystem.gass(37.70,1.570,44.2,0.781))
	print(CrudeOilSystem.gass(46.79,1.385,27.3,0.848))