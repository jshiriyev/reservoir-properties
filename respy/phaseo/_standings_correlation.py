import numpy as np

from _crude_oil_system import CrudeOilSystem as cos

class Standing:

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

		Cpb = p/18.2+1.4

		Rs = gg*(Cpb*10**x)**1.20482

		return Rs

	@staticmethod
	def bpp(T:float,API:float,gg:float,Rsb:float):
		"""
		Calculates the bubblepoint pressure of the oil at reservoir conditions
		
		Inputs:
		------
		T 	: System temperature, °F
		API : API gravity of oil, dimensionless
		gg	: Specific gravity of the separator gas
		Rsb	: Gas solubility at the bubble-point pressure, scf/STB

		The equations are valid to 325F. The correlation should be used with caution
		if nonhydrocarbon components are known to be present in the system.

		"""
		x = 0.00091*T-0.0125*API

		Cpb = (Rsb/gg)**0.83*10**x

		return 18.2*(Cpb-1.4)

	@staticmethod
	def fvf(T:float,API:float,gg:float,Rs:np.ndarray):
		"""
		Standing (1947) presented a graphical correlation for estimating the oil
		formation volume factor with the gas solubility, gas gravity, oil gravity,
		and reservoir temperature as the correlating parameters.

		This graphical correlation originated from examining a total of 105 experimental data
		points on 22 different California hydrocarbon systems. An average error
		of 1.2% was reported for the correlation.
		
		Standing (1981) showed that the oil formation volume factor can be expressed
		more conveniently in a mathematical form which requires following inputs:

		Inputs:
		------
		T 	: Temperature, °T
		API : API gravity of oil, dimensionless
		gg 	: Specific gravity of the solution gas
		Rs 	: Solution GOR, scf/STB

		"""
		go = cos.API2spgr(API)

		CBob = Rs*(gg/go)**0.5+1.25*T

		# oil formation volume factor at the bubble-point pressure, bbl/STB
		Bob = 0.9759+0.00012*CBob**1.2

        return Bob # Bob*np.exp(-co*(p-bpp))

	@staticmethod
	def comp(T:float,p:np.ndarray,API:float,gg:float,Rs:np.ndarray,fvfo:np.ndarray,fvfg:np.ndarray):
		"""
		It calculates compressibility based on the analytical derivation.

		Inputs:
		------
		T 	 : Temperature, °F
		p 	 : Pressure, psia
		API	 : API gravity of oil, dimensionless
		gg 	 : Specific gravity of the solution gas
		Rs 	 : Gas solubility at pressure p, scf/STB

		fvfo : Oil formation volume factor at p, bbl/STB
		fvfg : Gas formation volume factor at pressure p, bbl/scf

		"""
		go = cos.API2spgr(API)

		pRs = Rs/(0.83*p+21.75)
		pBo = 0.000144*pRs*np.sqrt(gg/go)*(Rs*np.sqrt(gg/go)+1.25*T)**0.12

		return (fvfg*pRs-pBo)/fvfo

if __name__ == "__main__":

	print(standing(250,2377+14.7,47.1,0.851))
	print(standing(220,2620+14.7,40.7,0.855))
	print(standing(260,2051+14.7,48.6,0.911))
	print(standing(237,2884+14.7,40.5,0.898))
	print(standing(218,3045+14.7,44.2,0.781))
	print(standing(180,4239+14.7,27.3,0.848))