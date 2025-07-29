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

if __name__ == "__main__":

	sgco = CrudeOilSystem.get_sgco(53)
	print(sgco)

	print(CrudeOilSystem.sgco_to_gAPI(sgco))

	print(CrudeOilSystem.get_sgsg((724,0.743),(202,0.956),ST=(58,1.296)))