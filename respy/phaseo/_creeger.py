class Creeger():

	def Rst():
		"""Calculates the stock-tank GOR. It should not be used
		if the separator temperature is >140F."""

		A1 =  0.3818
		A2 = -5.506
		A3 =  2.902
		A4 =  1.327
		A5 = -0.7355

		logrst = A1+A2*np.log(gamma_oil)+A3*np.log(gamma_gas_sp)+\
		            A4*np.log(Psp)+A5*np.log(Tsp)

		return np.exp(logrst)