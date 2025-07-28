import numpu as np

from ._crude_oil_system import CrudeOilSystem as cos

class PetroskyFarshad:

	@staticmethod
	def gass(T,p,API,gg)
		"""
		Petrosky and Farshad (1993) used a nonlinear multiple regression software
		to develop a gas solubility correlation. The authors constructed a
		PVT database from 81 laboratory analyses from the Gulf of Mexico crude
		oil system. Petrosky and Farshad proposed the following expression
		
		T = temperature, °F
		p = pressure, psia

		"""
		x = 7.916e-4*API**1.5410 - 4.561e-5*T**1.3911
		# print(f"\n{x=}")

		return ((p/112.27+12.340)*gg**0.8439*10**x)**1.73184

if __name__ == "__main__":

	print(petrosky_farshad(250,2377+14.7,47.1,0.851))
	print(petrosky_farshad(220,2620+14.7,40.7,0.855))
	print(petrosky_farshad(260,2051+14.7,48.6,0.911))
	print(petrosky_farshad(237,2884+14.7,40.5,0.898))
	print(petrosky_farshad(218,3045+14.7,44.2,0.781))
	print(petrosky_farshad(180,4239+14.7,27.3,0.848))