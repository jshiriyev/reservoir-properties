class GasOilInterfacialTension:
	
	# 25. surface/interfacial tension
	def oil_tens(P, T, oil_grav):
		"""Function to Calculate Gas-Oil Interfacial Tension in dynes/cm"""
		#P          pressure, psia
		#T          temperature, °F
		#oil_grav   API oil gravity
		s68 = 39 - 0.2571 * oil_grav
		s100 = 37.5 - 0.2571 * oil_grav
		if (T <= 68):
		    st = s68
		elif (T >= 100):
		    st = s100
		else:
		    st = s68 - (T - 68) * (s68 - s100) / 32

		c = 1 - 0.024 * P ** 0.45
		so = c * st
		if so < 1:
		    so = 1

		return so