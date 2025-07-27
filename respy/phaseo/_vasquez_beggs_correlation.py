import numpy as np

class VasquezBeggs:

	@staticmethod
	def gass(T,p,API,gg,Tsep,psep):
		"""
		Vasquez and Beggs (1980) presented an improved empirical correlation
		for estimating Rs. The correlation was obtained by regression analysis
		using 5,008 measured gas solubility data points.

		Based on oil gravity, the measured data were divided into two groups.
		This division was made at a value of oil gravity of 30°API.

		T		  temperature, °F
		P		  pressure, psia
		Tsep	   separator temperature, °F
		Psep	   separator pressure, psia
		Pb		 bubble point pressure, psia
		gas_grav   gas specific gravity
		oil_grav   API oil gravity

		An independent evaluation of the above correlation by Sutton and
		Farashad (1984) shows that the correlation is capable of predicting gas
		solubilities with an average absolute error of 12.7%.

		"""
		if API<30:
			C1 = 0.0362
			C2 = 1.0937
			C3 = 25.7240
		else:
			C1 = 0.0178
			C2 = 1.1870
			C3 = 23.931

		ggs = gg*(1+5.912e-5*API*Tsep*np.log10(psep/114.7))
		# print(f"\n{ggs=}")

		# gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)

		# if (P <= Pb):
		#	 Rs = C1 * gas_grav_corr * P** C2 * math.exp(C3 * oil_grav / (T + 460))
		# else:
		#	 Rs = C1 * gas_grav_corr * Pb ** C2 * math.exp(C3 * oil_grav / (T + 460))

		return C1*ggs*(p+14.7)**C2*np.exp(C3*API/(T+460))

	@staticmethod
	def bpp(T, Tsep, Psep, gas_grav, oil_grav, GOR):
		""" CFunction to Calculate Bubble Point Pressure in psia using Standing Correlation"""
		#T		  temperature, °F
		#Tsep	   separator temperature, °F
		#Psep	   separator pressure, psia
		#gas_grav   gas specific gravity
		#oil_grav   API oil gravity
		#GOR		producing gas-oil ratio, scf/stb
		gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
		if (oil_grav<= 30) :
			C1 = 0.0362
			C2 = 1.0937
			C3 = 25.724
		else:
			C1 = 0.0178
			C2 = 1.187
			C3 = 23.931
		
		Pbubl = (GOR / (C1 * gas_grav_corr * math.exp(C3 * oil_grav / (T + 460)))) **(1 / C2)
		return Pbubl

    # def get_isothermal_compressibility(self):
    # 18. isothermal compressibility coefficient of crude oil
    def comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav):
        """Function to Calculate Oil Isothermal Compressibility in 1/psi"""
        #'T          temperature, °F
        #'P          pressure, psia
        #'Tsep       separator temperature, °F
        #'Psep       separator pressure, psia
        #'Rs         solution gas-oil ratio, scf/stb
        #'gas_grav   gas specific gravity
        #'oil_grav   API oil gravity

        # A1 =  -1_433.0
        # A2 =       5.0
        # A3 =      17.2
        # A4 =  -1_180.0
        # A5 =      12.61
        # A6 = 100_000.0

        # return (A1+A2*Rs+A3*temp+A4*gamma_gas+A5*gamma_API)/(A6*pres)
        
        gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
        oil_compr = (5 * Rs + 17.2 * T - 1180 * gas_grav_corr + 12.61 * oil_grav - 1433) / (P * 10 ** 5)
        return oil_compr

   	# def get_fvf(self):
	# 17. oil formation volume factor
	# 19. oil formation volume factor for undersaturated oils
	def fvf(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
		"""Function to Calculate Oil Formation Volume Factor in bbl/stb"""
		#'T		  temperature, °F
		#P		  pressure, psia
		#Tsep	   separator temperature, °F
		#Psep	   separator pressure, psia
		#Pb		 bubble point pressure, psia
		#Rs		 solution gas-oil ratio, scf/stb
		#gas_grav   gas specific gravity
		#oil_grav   API oil gravity
		gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
		if (oil_grav <= 30) :
			C1 = 0.0004677
			C2 = 1.751E-05
			C3 = -1.811E-08
		else:
			C1 = 0.000467
			C2 = 1.1E-05
			C3 = 1.337E-09
		
		if (P <= Pb):
			Bo = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr) + C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
		else:
			Bob = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr)+ C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
			co = oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav)
			Bo = Bob * math.exp(co * (Pb - P))
		
		return  Bo
		
if __name__ == "__main__":

	print(vasquez_beggs(250,2377,47.1,0.851, 60,150+14.7))
	print(vasquez_beggs(220,2620,40.7,0.855, 75,100+14.7))
	print(vasquez_beggs(260,2051,48.6,0.911, 72,100+14.7))
	print(vasquez_beggs(237,2884,40.5,0.898,120, 60+14.7))
	print(vasquez_beggs(218,3045,44.2,0.781, 60,200+14.7))
	print(vasquez_beggs(180,4239,27.3,0.848,173, 85+14.7))