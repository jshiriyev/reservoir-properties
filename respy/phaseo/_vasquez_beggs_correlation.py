import numpy as np

from ._crude_oil_system import VasquezBeggs as cos

class VasquezBeggs:

	@staticmethod
	def bpp(T:float,API:float,sgsg:float,Rsb:float,Tsep:float=None,psep:float=None):
		"""
		Vasquez and Beggs’ gas solubility correlation can be solved for the
		bubble-point pressure that will require the following inputs:
		
		Inputs:
		------
		T	 : System temperature, °F
		API  : API oil gravity, dimensionless
		sgsg : Specific gravity of the separator gas
		Rsb	 : Gas solubility at the bubble-point pressure, scf/STB

		Tsep : Separator temperature, °F
		psep : Separator pressure, psia

		"""
		if API<=30.:
			C1 = 27.624
			C2 = 0.914328
			C3 = 11.172
		else:
			C1 = 56.18
			C2 = 0.84246
			C3 = 10.393

		if Tsep and psep:
			sgsg = VasquezBeggs.sgsgcorr(Tsep,psep,API,sgsg)
		
		a = -C3*API/T

		# return (Rs / (C1 * gg * math.exp(C3 * API / (T + 460)))) **(1 / C2)
		return (C1*Rsb/sgsg*10**a)**C2

	@staticmethod
	def gass(T:float,p:np.ndarray,API:float,sgsg:float,Tsep:float=None,psep:float=None):
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
		sgsg 		: Solution gas specific gravity

		Tsep	: Separator temperature, °F
		psep	: Separator pressure, psia

		An independent evaluation of the above correlation by Sutton and
		Farashad (1984) shows that the correlation is capable of predicting gas
		solubilities with an average absolute error of 12.7%.

		"""
		if API<=30.:
			C1 = 0.0362
			C2 = 1.0937
			C3 = 25.7240
		else:
			C1 = 0.0178
			C2 = 1.1870
			C3 = 23.931

		if Tsep and psep:
			sgsg = VasquezBeggs.sgsgcorr(Tsep,psep,API,sgsg)

		return C1*sgsg*p**C2*np.exp(C3*API/(T+460))

	def fvf(T:float,API:float,sgsg:float,Rs:np.ndarray,Tsep:float=None,psep:float=None):
		"""Calculates Oil Formation Volume Factor in bbl/stb

		Vasquez and Beggs (1980) developed a relationship for determining
		Bo as a function of Rs, go, sgsg, and T. The proposed correlation was based
		on 6,000 measurements of Bo at various pressures. Using the regression
		analysis technique, Vasquez and Beggs found the method to
		be the best form to reproduce the measured data
		
		Inputs:
		------
		T	 : temperature, °F
		API  : API oil gravity
		sgsg : gas specific gravity
		Rs	 : solution gas-oil ratio, scf/stb

		Tsep : separator temperature, °F
		psep : separator pressure, psia

		Vasquez and Beggs reported an average error of 4.7% for the proposed
		correlation.

		"""
		if API<=30.:
			C1 = 4.677E-04
			C2 = 1.751E-05
			C3 = -1.811E-08
		else:
			C1 = 4.670E-4
			C2 = 1.100E-05
			C3 = 1.337E-09

		if Tsep and psep:
			sgsg = VasquezBeggs.sgsgcorr(Tsep,psep,API,sgsg)
		
		return 1.+C1*Rs+(T-60)*(API/sgsg)*(C2+C3*Rs)
		
	def comp(T:float,p:np.ndarray,API:float,sgsg:float,Rsb:float,Tsep:float=None,psep:float=None):
		"""Calculates oil isothermal compressibility in 1/psi

		From a total of 4,036 experimental data points used in a linear regression
		model, Vasquez and Beggs (1980) correlated the isothermal oil compressibility
		coefficients with Rsb, T, °API, gg, and p. They proposed the
		following expression:
		
		Inputs:
		------
		T 	 : temperature, °F
		p  	 : pressure, psia
		API  : API oil gravity
		sgsg : gas specific gravity
		Rsb  : solution gas-oil ratio, scf/stb

		Tsep : separator temperature, °F
		psep : separator pressure, psia

		"""
		if Tsep and psep:
			sgsg = VasquezBeggs.sgsgcorr(Tsep,psep,API,sgsg)

		return (-1433.+5.*Rsb+17.2*T-1180.*sgsg+12.61*API)/(p*10**5)

	@staticmethod
	def sgsgcorr(Tsep,psep,API,sgsg):
		"""Method to calculate corrected specific gravity for the solution gas:

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
		sgsg  	: Gas gravity at the actual separator conditions of psep and Tsep

		Returns:
		-------
		Gas gravity at the reference separator pressure

		The gas gravity used to develop all the correlations reported by the
		authors was that which would result from a two-stage separation. The
		first-stage pressure was chosen as 100 psig and the second stage was the
		stock tank. If the separator conditions are unknown, the unadjusted gas
		gravity may be used.

		"""
		return sgsg*(1+5.912e-5*API*Tsep*np.log10(psep/(100+14.7)))
		
if __name__ == "__main__":

	print(VasquezBeggs.gass(250,2377,47.1,0.851, 60,150+14.7))
	print(VasquezBeggs.gass(220,2620,40.7,0.855, 75,100+14.7))
	print(VasquezBeggs.gass(260,2051,48.6,0.911, 72,100+14.7))
	print(VasquezBeggs.gass(237,2884,40.5,0.898,120, 60+14.7))
	print(VasquezBeggs.gass(218,3045,44.2,0.781, 60,200+14.7))
	print(VasquezBeggs.gass(180,4239,27.3,0.848,173, 85+14.7))

	print(VasquezBeggs.sgsgcorr( 60,150+14.7,47.1,0.851))
	print(VasquezBeggs.sgsgcorr( 75,100+14.7,40.7,0.855))
	print(VasquezBeggs.sgsgcorr( 72,100+14.7,48.6,0.911))
	print(VasquezBeggs.sgsgcorr(120, 60+14.7,40.5,0.898))
	print(VasquezBeggs.sgsgcorr( 60,200+14.7,44.2,0.781))
	print(VasquezBeggs.sgsgcorr(173, 85+14.7,27.3,0.848))