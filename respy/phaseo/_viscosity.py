class Viscosity:

	# def get_viscosity(self):
    # 21. crude oil viscosity
    # 22. methods of calculating the dead oil viscosity
    # 23. methods of calculating the saturated oil viscosity
    # 24. methods of calculating the viscosity of the undersaturated oil
    def oil_visc(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
        """Function to Calculate Oil Viscosity in cp"""
        #'T          temperature, °F
        #'P          pressure, psia
        #'Tsep       separator temperature, °F
        #'Psep       separator pressure, psia
        #'Pb         bubble point pressure, psia
        #'Rs         solution gas-oil ratio, scf/stb
        #'gas_grav   gas specific gravity
        #'oil_grav   API oil gravity
        
        a = 10.715 * (Rs + 100) ** (-0.515)
        b = 5.44 * (Rs + 150) ** (-0.338)
        Z = 3.0324 - 0.0203 * oil_grav
        Y = 10 **Z
        x = Y * T ** (-1.163)
        visc_oD = 10 ** x - 1
        if (P <= Pb):
            visc_o = a * visc_oD ** b
        else:
            M = 2.6 * P ** 1.187 * math.exp(-11.513 - 8.98E-05 * P)
            visc_ob = a * visc_oD ** b
            visc_o = visc_ob * (P / Pb) ** M

        return visc_o