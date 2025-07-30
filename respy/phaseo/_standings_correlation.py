import numpy as np

from ._crude_oil_system import CrudeOilSystem as cos

class StandingsCorrelation:

	@staticmethod
	def bpp(Rsb:float,sgsg:float,gAPI:float,temp:float):
		"""
		Calculates the bubblepoint pressure of the oil at reservoir conditions
		
		Inputs:
		------
		Rsb	 : Gas solubility at the bubble-point pressure, scf/STB
		
		sgsg : Specific gravity of the separator gas
		gAPI : API gravity of oil, dimensionless
		temp : System temperature, °F

		The equations are valid to 325F. The correlation should be used with caution
		if nonhydrocarbon components are known to be present in the system.

		"""
		x = 0.00091*temp-0.0125*gAPI

		Cpb = (Rsb/sgsg)**0.83*10**x

		return 18.2*(Cpb-1.4)

	@staticmethod
	def gass(p:np.ndarray,bpp:float,sgsg:float,gAPI:float,temp:float):
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
		p 	 : System pressure, psia

		bpp  : Bubble point pressure, psia
		
		sgsg : Solution gas specific gravity
		gAPI : API gravity of oil, dimensionless
		temp : System temperature, °F

		It should be noted that Standing’s equation is valid for applications at
		and below the bubble-point pressure of the crude oil.

		"""
		Rsb = StandingsCorrelation.gass_sat(bpp,sgsg,gAPI,temp)

		_p = np.atleast_1d(p)
		Rs = np.full_like(_p,Rsb)

		Rs[_p<bpp] = StandingsCorrelation.gass_sat(_p[_p<bpp],sgsg,gAPI,temp)

		return Rs

	@staticmethod
	def gass_sat(p:float|np.ndarray,sgsg:float,gAPI:float,temp:float):

		x = 0.0125*gAPI-0.00091*temp

		return sgsg*((p/18.2+1.4)*10**x)**(1./0.83)

	@staticmethod
	def gass_sat_prime(p:float|np.ndarray,sgsg:float,gAPI:float,temp:float):

		sgco = cos.gAPI_to_sgco(gAPI)

		gass = StandingsCorrelation.gass_sat(p,sgsg,gAPI,temp)

		return gass/(0.83*p+21.75)

	@staticmethod
	def fvf_sat(p:float|np.ndarray,sgsg:float,gAPI:float,temp:float):
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
		gass : Solution GOR, scf/STB
		
		sgsg : Specific gravity of the solution gas
		gAPI : API gravity of oil, dimensionless
		temp : Temperature, °T

		"""
		sgco = cos.gAPI_to_sgco(gAPI)

		gass = StandingsCorrelation.gass_sat(p,sgsg,gAPI,temp)

		CBob = gass*(sgsg/sgco)**0.5+1.25*temp

		return 0.9759+0.00012*CBob**1.2

	@staticmethod
	def fvf_sat_prime(p:float|np.ndarray,sgsg:float,gAPI:float,temp:float):

		sgco = cos.gAPI_to_sgco(gAPI)
		gass = StandingsCorrelation.gass_sat(p,sgsg,gAPI,temp)

		return 0.000144*(gass/(0.83*p+21.75))*np.sqrt(sgsg/sgco)*(gass*np.sqrt(sgsg/sgco)+1.25*temp)**0.12

	@staticmethod
	def comp_sat(p:float|np.ndarray,fvfg:float|np.ndarray,fvfo:float|np.ndarray,sgsg:float,gAPI:float,temp:float):
		"""
		It calculates compressibility based on the analytical derivation.

		Inputs:
		------
		p 	 : Pressure, psia

		fvfg : Gas formation volume factor at pressure p, bbl/scf
		fvfo : Oil formation volume factor at p, bbl/STB

		sgsg : Specific gravity of the solution gas
		gAPI : API gravity of oil, dimensionless
		temp : Temperature, °F

		"""
		sgco = cos.gAPI_to_sgco(gAPI)

		gass = StandingsCorrelation.gass(p,bpp,sgsg,gAPI,temp)

		pRs = gass/(0.83*p+21.75)
		pBo = 0.000144*pRs*np.sqrt(sgsg/sgco)*(gass*np.sqrt(sgsg/sgco)+1.25*temp)**0.12

		return (fvfg*pRs-pBo)/fvfo

if __name__ == "__main__":

	print(StandingsCorrelation.gass(250,2377+14.7,47.1,0.851))
	print(StandingsCorrelation.gass(220,2620+14.7,40.7,0.855))
	print(StandingsCorrelation.gass(260,2051+14.7,48.6,0.911))
	print(StandingsCorrelation.gass(237,2884+14.7,40.5,0.898))
	print(StandingsCorrelation.gass(218,3045+14.7,44.2,0.781))
	print(StandingsCorrelation.gass(180,4239+14.7,27.3,0.848))