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
	def spgr(rhoo:float|np.ndarray,rhow:float=62.4):
		"""
		The specific gravity of a crude oil is defined as the ratio of the
		density of the oil to that of water. Both densities are measured
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
	def API(spgr:float|np.ndarray):
		"""
		Calculates API gravity from specific gravity:
		
		Inputs:
		------
		spgr	: specific gravity of the oil

		Returns:
		-------
		API gravity of the oil

		The API gravities of crude oils usually range from 47° API for the
		lighter crude oils to 10° API for the heavier asphaltic crude oils.

		"""
		return 141.5/spgr-131.5

	@staticmethod
	def spgr_solg(*separators,ST:tuple[float,float]):
		"""
		Calculates the specific gravity of the solution gas.

		The specific gravity of the solution gas is described by the weighted
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

	@staticmethod
	def gas_correction(Tsep,psep,API,gg):
		"""Method to calculate corrected gas gravity:

		Realizing that the value of the specific gravity of the gas depends on
		the conditions under which it is separated from the oil, Vasquez and
		Beggs proposed that the value of the gas specific gravity as obtained
		from a separator pressure of 100 psig be used in the above equation. This
		reference pressure was chosen because it represents the average field
		separator conditions.
		
		Inputs:
		------
		Tsep	: Actual separator temperature, °F
		psep	: Actual separator pressure, psia
		API	 	: API oil gravity
		gg	  	: Gas gravity at the actual separator conditions of psep and Tsep

		Returns:
		-------
		Gas gravity at the reference separator pressure

		The gas gravity used to develop all the correlations reported by the
		authors was that which would result from a two-stage separation. The
		first-stage pressure was chosen as 100 psig and the second stage was the
		stock tank. If the separator conditions are unknown, the unadjusted gas
		gravity may be used.

		"""
		return gg*(1+5.912e-5*API*Tsep*np.log10(psep/(100+14.7)))

if __name__ == "__main__":

	sg = CrudeOilSystem.spgr(53)
	print(sg)

	print(CrudeOilSystem.API(sg))

	print(CrudeOilSystem.spgr_solg((724,0.743),(202,0.956),ST=(58,1.296)))

	print(CrudeOilSystem.gas_correction( 60,150+14.7,47.1,0.851))
	print(CrudeOilSystem.gas_correction( 75,100+14.7,40.7,0.855))
	print(CrudeOilSystem.gas_correction( 72,100+14.7,48.6,0.911))
	print(CrudeOilSystem.gas_correction(120, 60+14.7,40.5,0.898))
	print(CrudeOilSystem.gas_correction( 60,200+14.7,44.2,0.781))
	print(CrudeOilSystem.gas_correction(173, 85+14.7,27.3,0.848))