import numpy as np

from _crude_oil_system import CrudeOilSystem as cos

class VasquezBeggs:

	@staticmethod
	def gass(T:float,p:np.ndarray,API:float,gg:float,Tsep:float=None,psep:float=None):
		"""
		Vasquez and Beggs (1980) presented an improved empirical correlation
		for estimating Rs. The correlation was obtained by regression analysis
		using 5,008 measured gas solubility data points.

		Based on oil gravity, the measured data were divided into two groups.
		This division was made at a value of oil gravity of 30°API.
		
		Inputs:
		------
		T		: System temperature, °F
		p		: System pressure, psia
		API 	: API gravity of oil, dimensionless
		gg 		: Solution gas specific gravity

		Tsep	: Separator temperature, °F
		psep	: Separator pressure, psia

		An independent evaluation of the above correlation by Sutton and
		Farashad (1984) shows that the correlation is capable of predicting gas
		solubilities with an average absolute error of 12.7%.

		"""
		if API<=30:
			C1 = 0.0362
			C2 = 1.0937
			C3 = 25.7240
		else:
			C1 = 0.0178
			C2 = 1.1870
			C3 = 23.931

		if Tsep and psep:
			gg = cos.gas_correction(Tsep,psep,API,gg)

		return C1*gg*p**C2*np.exp(C3*API/(T+460))

	@staticmethod
	def bpp(T:float,API:float,gg:float,Rsb:float,Tsep:float=None,psep:float=None):
		"""
		Vasquez and Beggs’ gas solubility correlation can be solved for the
		bubble-point pressure that will require the following inputs:
		
		Inputs:
		------
		T	 : System temperature, °F
		API  : API oil gravity, dimensionless
		gg	 : Specific gravity of the separator gas
		Rsb	 : Gas solubility at the bubble-point pressure, scf/STB

		Tsep : Separator temperature, °F
		psep : Separator pressure, psia

		"""
		if API<=30:
			C1 = 27.624
			C2 = 0.914328
			C3 = 11.172
		else:
			C1 = 56.18
			C2 = 0.84246
			C3 = 10.393

		if Tsep and psep:
			gg = cos.gas_correction(Tsep,psep,API,gg)
		
		a = -C3*API/T

		# return (Rs / (C1 * gg * math.exp(C3 * API / (T + 460)))) **(1 / C2)
		return (C1*Rsb/gg*10**a)**C2

	def fvf(T,p,Tsep,psep,Pb,Rs,gg,API):
		"""Function to Calculate Oil Formation Volume Factor in bbl/stb
		
		Inputs:
		------
		T		temperature, °F
		p		pressure, psia
		Tsep	separator temperature, °F
		psep	separator pressure, psia
		Pb		bubble point pressure, psia
		Rs		solution gas-oil ratio, scf/stb
		API 	API oil gravity
		gg   	gas specific gravity

		"""
		gg_corr = correct(Tsep, psep, gg, API)

		if (API <= 30) :
			C1 = 0.0004677
			C2 = 1.751E-05
			C3 = -1.811E-08
		else:
			C1 = 0.000467
			C2 = 1.1E-05
			C3 = 1.337E-09
		
		if (P <= Pb):
			Bo = 1 + C1 * Rs + C2 * (T - 60) * (API / gg_corr) + C3 * Rs * (T - 60) * (API / gg_corr)
		else:
			Bob = 1 + C1 * Rs + C2 * (T - 60) * (API / gg_corr)+ C3 * Rs * (T - 60) * (API / gg_corr)
			co = oil_comp(T, P, Tsep, psep, Rs, gg, API)
			Bo = Bob * math.exp(co * (Pb - P))
		
		return  Bo

	def comp(T,p,Tsep,psep,Rs,gg,API):
		"""Function to Calculate Oil Isothermal Compressibility in 1/psi
		
		Inputs:
		------
		T 		temperature, °F
		p 		pressure, psia
		Tsep	separator temperature, °F
		psep	separator pressure, psia
		Rs 		solution gas-oil ratio, scf/stb
		API   	API oil gravity
		gg   	gas specific gravity

		"""
		# A1 =  -1_433.0
		# A2 =	   5.0
		# A3 =	  17.2
		# A4 =  -1_180.0
		# A5 =	  12.61
		# A6 = 100_000.0

		# return (A1+A2*Rs+A3*temp+A4*gamma_gas+A5*gamma_API)/(A6*pres)
		
		gg_corr = correct(Tsep, psep, gg, API)

		oil_compr = (5 * Rs + 17.2 * T - 1180 * gg_corr + 12.61 * API - 1433) / (P * 10 ** 5)

		return oil_compr
		
if __name__ == "__main__":

	print(vasquez_beggs(250,2377,47.1,0.851, 60,150+14.7))
	print(vasquez_beggs(220,2620,40.7,0.855, 75,100+14.7))
	print(vasquez_beggs(260,2051,48.6,0.911, 72,100+14.7))
	print(vasquez_beggs(237,2884,40.5,0.898,120, 60+14.7))
	print(vasquez_beggs(218,3045,44.2,0.781, 60,200+14.7))
	print(vasquez_beggs(180,4239,27.3,0.848,173, 85+14.7))