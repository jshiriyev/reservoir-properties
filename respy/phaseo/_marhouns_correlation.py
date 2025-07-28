import numpy as np

from ._crude_oil_system import CrudeOilSystem as cos

class Marhoun:

	@staticmethod
	def gass(p,sgsg,gAPI,temp):
		"""
		Marhoun (1988) developed an expression for estimating the saturation
		pressure of the Middle Eastern crude oil systems. The correlation originates
		from 160 experimental saturation pressure data.
		
		p 	 = system pressure, psia
		
		sgsg = gas specific gravity
		sgco = stock-tank oil gravity
		temp = temperature, °F
		
		"""
		a = 185.843208
		b = 1.877840
		c = -3.1437
		d = -1.32657
		e = 1.398441

		sgco = 141.5/(gAPI+131.5)

		return (a*sgsg**b*sgco**c*(temp+460)**d*p)**e

if __name__ == "__main__":

	print(marhoun(250,2377+14.7,47.1,0.851))
	print(marhoun(220,2620+14.7,40.7,0.855))
	print(marhoun(260,2051+14.7,48.6,0.911))
	print(marhoun(237,2884+14.7,40.5,0.898))
	print(marhoun(218,3045+14.7,44.2,0.781))
	print(marhoun(180,4239+14.7,27.3,0.848))