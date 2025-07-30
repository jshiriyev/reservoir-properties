import numpy as np

from ._crude_oil_system import CrudeOilSystem as cos

class VasquezBeggsCorrelation:

	@staticmethod
	def bpp(Rsb:float,sgsg:float,gAPI:float,temp:float,psep:float=None,Tsep:float=None):
		"""
		Vasquez and Beggs’ gas solubility correlation can be solved for the
		bubble-point pressure that will require the following inputs:
		
		Inputs:
		------
		Rsb	 : Gas solubility at the bubble-point pressure, scf/STB
		
		sgsg : Specific gravity of the separator gas
		gAPI : API oil gravity, dimensionless
		temp : System temperature, °F
		
		psep : Separator pressure, psia
		Tsep : Separator temperature, °F

		"""
		if gAPI<=30.:
			C1 = 27.624
			C2 = 0.914328
			C3 = 11.172
		else:
			C1 = 56.18
			C2 = 0.84246
			C3 = 10.393

		sgsg = VasquezBeggsCorrelation.sgsgcorr(sgsg,gAPI,psep,Tsep)
		
		a = -C3*gAPI/temp

		# return (Rs / (C1 * sgsg * math.exp(C3 * gAPI / (temp + 460)))) **(1 / C2)
		return (C1*Rsb/sgsg*10**a)**C2

	@staticmethod
	def gass(p:np.ndarray,bpp:float,sgsg:float,gAPI:float,temp:float,psep:float=None,Tsep:float=None):
		"""
		Vasquez and Beggs (1980) presented an improved empirical correlation
		for estimating Rs. The correlation was obtained by regression analysis
		using 5,008 measured gas solubility data points.

		Based on oil gravity, the measured data were divided into two groups.
		This division was made at a value of oil gravity of 30°API.
		
		Inputs:
		------
		p		: System pressure, psia

		bpp		: Bubble point pressure, psia
		
		sgsg	: Solution gas specific gravity
		gAPI 	: API gravity of oil, dimensionless
		temp	: System temperature, °F
		
		psep	: Separator pressure, psia
		Tsep	: Separator temperature, °F

		An independent evaluation of the above correlation by Sutton and
		Farashad (1984) shows that the correlation is capable of predicting gas
		solubilities with an average absolute error of 12.7%.

		"""
		if gAPI<=30.:
			C1 = 0.0362
			C2 = 1.0937
			C3 = 25.7240
		else:
			C1 = 0.0178
			C2 = 1.1870
			C3 = 23.931

		sgsg = VasquezBeggsCorrelation.sgsgcorr(sgsg,gAPI,psep,Tsep)

		p = np.atleast_1d(p)

		Rsb = C1*sgsg*bpp**C2*np.exp(C3*gAPI/(temp+460))

		_Rs = np.full_like(p,Rsb)
		_Rs[p<bpp] = C1*sgsg*p[p<bpp]**C2*np.exp(C3*gAPI/(temp+460))

		return _Rs

	@staticmethod
	def fvf(p:np.ndarray,bpp:float,sgsg:float,gAPI:float,temp:float,psep:float=None,Tsep:float=None):
		"""Calculates Oil Formation Volume Factor in bbl/stb

		Vasquez and Beggs (1980) developed a relationship for determining
		Bo as a function of Rs, go, sgsg, and T. The proposed correlation was based
		on 6,000 measurements of Bo at various pressures. Using the regression
		analysis technique, Vasquez and Beggs found the method to
		be the best form to reproduce the measured data
		
		Inputs:
		------
		p	 : System pressure, psia

		bpp	 : Bubble point pressure, psia
		
		sgsg : gas specific gravity
		gAPI : API oil gravity
		temp : temperature, °F
		
		psep : separator pressure, psia
		Tsep : separator temperature, °F

		Vasquez and Beggs reported an average error of 4.7% for the proposed
		correlation.

		"""
		if gAPI<=30.:
			C1 = 4.677E-04
			C2 = 1.751E-05
			C3 = -1.811E-08
		else:
			C1 = 4.670E-4
			C2 = 1.100E-05
			C3 = 1.337E-09

		sgsg = VasquezBeggsCorrelation.sgsgcorr(sgsg,gAPI,psep,Tsep)

		p = np.atleast_1d(p)

		Rsb = VasquezBeggsCorrelation.gass(bpp,bpp,sgsg,gAPI,temp)
		Bob = 1.+C1*Rsb+(temp-60)*(gAPI/sgsg)*(C2+C3*Rsb)

		Rs = VasquezBeggsCorrelation.gass(p,bpp,sgsg,gAPI,temp)
		Bo = 1.+C1*Rs+(temp-60)*(gAPI/sgsg)*(C2+C3*Rs)

		A = (-1433.+5.*Rsb+17.2*temp-1180.*sgsg+12.61*gAPI)/(10**5)

		Bo[p>bpp] = Bob*np.exp(-A*np.log(p[p>bpp]/bpp))
		
		return Bo
	
	@staticmethod
	def comp(p:np.ndarray,bpp:float,fvfg:np.ndarray,fvfo:np.ndarray,sgsg:float,gAPI:float,temp:float,psep:float=None,Tsep:float=None):
		"""Calculates oil isothermal compressibility in 1/psi

		From a total of 4,036 experimental data points used in a linear regression
		model, Vasquez and Beggs (1980) correlated the isothermal oil compressibility
		coefficients with Rsb, T, °API, gg, and p. They proposed the
		following expression:
		
		Inputs:
		------
		p  	 : pressure, psia

		Rsb  : solution gas-oil ratio, scf/stb
		
		sgsg : gas specific gravity
		gAPI : API oil gravity
		temp : temperature, °F
		
		psep : separator pressure, psia
		Tsep : separator temperature, °F

		"""

		# IT IS ONLY VALID FOR PRESSURES ABOVE BPP
		
		sgsg = VasquezBeggsCorrelation.sgsgcorr(sgsg,gAPI,psep,Tsep)

		Rsb = VasquezBeggsCorrelation.gass(bpp,bpp,sgsg,gAPI,temp)

		return (-1433.+5.*Rsb+17.2*temp-1180.*sgsg+12.61*gAPI)/(p*10**5)

	@staticmethod
	def sgsgcorr(sgsg:float,gAPI:float,psep:float=None,Tsep:float=None):
		"""Method to calculate corrected specific gravity for the solution gas:

		Realizing that the value of the specific gravity of the gas depends on
		the conditions under which it is separated from the oil, Vasquez and
		Beggs proposed that the value of the gas specific gravity as obtained
		from a separator pressure of 100 psig be used in the above equation. This
		reference pressure was chosen because it represents the average field
		separator conditions.
		
		Inputs:
		------
		sgsg  	: Gas gravity at the actual separator conditions of psep and Tsep
		gAPI	: API oil gravity
		psep	: Actual separator pressure, psia
		Tsep	: Actual separator temperature, °F

		Returns:
		-------
		Gas gravity at the reference separator pressure

		The gas gravity used to develop all the correlations reported by the
		authors was that which would result from a two-stage separation. The
		first-stage pressure was chosen as 100 psig and the second stage was the
		stock tank. If the separator conditions are unknown, the unadjusted gas
		gravity may be used.

		"""
		if psep and Tsep:
			return sgsg*(1+5.912e-5*gAPI*Tsep*np.log10(psep/(100+14.7)))

		return sgsg
		
if __name__ == "__main__":

	print(VasquezBeggsCorrelation.gass(250,2377,47.1,0.851, 60,150+14.7))
	print(VasquezBeggsCorrelation.gass(220,2620,40.7,0.855, 75,100+14.7))
	print(VasquezBeggsCorrelation.gass(260,2051,48.6,0.911, 72,100+14.7))
	print(VasquezBeggsCorrelation.gass(237,2884,40.5,0.898,120, 60+14.7))
	print(VasquezBeggsCorrelation.gass(218,3045,44.2,0.781, 60,200+14.7))
	print(VasquezBeggsCorrelation.gass(180,4239,27.3,0.848,173, 85+14.7))

	print(VasquezBeggsCorrelation.sgsgcorr( 60,150+14.7,47.1,0.851))
	print(VasquezBeggsCorrelation.sgsgcorr( 75,100+14.7,40.7,0.855))
	print(VasquezBeggsCorrelation.sgsgcorr( 72,100+14.7,48.6,0.911))
	print(VasquezBeggsCorrelation.sgsgcorr(120, 60+14.7,40.5,0.898))
	print(VasquezBeggsCorrelation.sgsgcorr( 60,200+14.7,44.2,0.781))
	print(VasquezBeggsCorrelation.sgsgcorr(173, 85+14.7,27.3,0.848))