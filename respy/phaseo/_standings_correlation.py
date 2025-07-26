import numpy as np

from ._crude_oil_system import CrudeOilSystem

class Standing(CrudeOilSystem):

	@staticmethod
	def gass(T:float,p:np.ndarray,API:float,gg:float):
		"""
		Standing (1947) proposed a graphical correlation for determining the
		gas solubility as a function of pressure, gas specific gravity, API gravity,
		and system temperature.

		The correlation was developed from a total of 105 experimentally determined
		data points on 22 hydrocarbon mixtures from California crude oils and
		natural gases. The proposed correlation has an average error of 4.8%.
		
		Standing (1981) expressed his proposed graphical correlation in the
		more convenient mathematical form which requires following inputs:
		
		Inputs:
		------
		T 	: System temperature, °F
		p 	: System pressure, psia
		API : API gravity of oil, dimensionless
		gg 	: Solution gas specific gravity

		It should be noted that Standing’s equation is valid for applications at
		and below the bubble-point pressure of the crude oil.

		"""
		x = 0.0125*API-0.00091*T
		# print(f"\n{x=}")

		Cpb = p/18.2+1.4

		Rs = gg*(Cpb*10**x)**1.20482

		return Rs

	@staticmethod
	def bpp(T,API,gg,Rs):
		"""
		Calculates the bubblepoint pressure of the oil at reservoir conditions
		
		T 	: System temperature, °F
		API : API gravity of oil, dimensionless
		gg	: specific gravity of the separator gas
		Rs	: solution GOR including stock-tank vent gas

		The equations are valid to 325F. The correlation should be used with caution
		if nonhydrocarbon components are known to be present in the system.

		"""
		x = 0.00091*T-0.0125*API

		Cpb = (Rs/gg)**0.83*10**x

		return 18.2*(Cpb-1.4)

	@staticmethod
	def fvf(T,API,gg,Rs):
		"""
		Standing (1947) presented a graphical correlation for estimating the oil
		formation volume factor with the gas solubility, gas gravity, oil gravity,
		and reservoir temperature as the correlating parameters.

		This graphical correlation originated from examining a total of 105 experimental data
		points on 22 different California hydrocarbon systems. An average error
		of 1.2% was reported for the correlation.
		
		Standing (1981) showed that the oil formation volume factor can be
		expressed more conveniently in a mathematical form by the following
		equation:

		"""
		# def fvf(pres,temp,Rs,gamma_gas,gamma_oil):

        # CBob = Rs*(gamma_gas/gamma_oil)**0.5+1.25*temp

        # Bob = 0.9759+12e-5*CBob**1.2

        # return Bob*numpy.exp(co*(pbubble-pres))
        
		return 0.9759+0.00012*(Rs*(gg/go)**0.5+1.25*T)**1.2

	@staticmethod
	def comp():
		"""


		"""
		pass

	@staticmethod
	def rho():
		"""


		"""
		pass

if __name__ == "__main__":

	print(standing(250,2377+14.7,47.1,0.851))
	print(standing(220,2620+14.7,40.7,0.855))
	print(standing(260,2051+14.7,48.6,0.911))
	print(standing(237,2884+14.7,40.5,0.898))
	print(standing(218,3045+14.7,44.2,0.781))
	print(standing(180,4239+14.7,27.3,0.848))