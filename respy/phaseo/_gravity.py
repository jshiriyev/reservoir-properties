"""
Crude Oil Gravity

The crude oil density is defined as the mass of a unit volume of the
crude at a specified pressure and temperature. It is usually expressed in
pounds per cubic foot. 

"""

import numpy as np

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

def API(spgr:float|np.ndarray):
	"""
	Calculates API gravity from specific gravity:
	
	Inputs:
	------
	spgr 	: specific gravity of the oil

	Returns:
	-------
	API gravity of the oil

	The API gravities of crude oils usually range from 47° API for the
	lighter crude oils to 10° API for the heavier asphaltic crude oils.

	"""
	return 141.5/spgr-131.5

def solg(stock_tank:tuple[float,float],*separators):
	"""
	Calculates the specific gravity of the solution gas.

	The specific gravity of the solution gas is described by the weighted
	average of the specific gravities of the separated gas from each separator.

	This weighted-average approach is based on the separator gas-oil ratio.
	
	The stock_tank is a tuple with Rst and gst values:
		Rst = gas-oil ratio from the stock tank, scf/ STB
		gst = gas gravity from the stock tank
	
	For each separator Rsep and gsep values must be provided:
		Rsep = separator gas-oil ratio, scf/STB
		gsep = separator gas gravity
	
	Returns:
	-------
	Specific gravity of the solution gas.

	"""
	Rst,gst = stock_tank

	upper = Rst*gst
	lower = Rst

	for sep in separators:
		Rsep,gsep = sep

		upper += Rsep*gsep
		lower += Rsep

	return upper/lower

if __name__ == "__main__":

	sg = spgr(53)

	print(sg)
	print(API(sg))

	print(solg((58,1.296),(724,0.743),(202,0.956)))