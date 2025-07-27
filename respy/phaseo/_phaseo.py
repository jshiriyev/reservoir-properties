import os, sys

import numpy as np

from ._crude_oil_system import CrudeOilSystem

class phaseo(CrudeOilSystem):

    def gass(self,):
        """
        The gas solubility Rs is defined as the number of standard cubic feet of
        gas which will dissolve in one stock-tank barrel of crude oil at certain
        pressure and temperature. The solubility of a natural gas in a crude oil is a
        strong function of the pressure, temperature, API gravity, and gas gravity.

        For a particular gas and crude oil to exist at a constant temperature, the
        solubility increases with pressure until the saturation pressure is reached.
        At the saturation pressure (bubble-point pressure) all the available gases
        are dissolved in the oil and the gas solubility reaches its maximum value.
        
        Rather than measuring the amount of gas that will dissolve in a given
        stock-tank crude oil as the pressure is increased, it is customary to determine
        the amount of gas that will come out of a sample of reservoir crude
        oil as pressure decreases.

        A typical gas solubility curve, as a function of pressure for an undersaturated
        crude oil:

        - As the pressure is reduced from the initial reservoir pressure pi, to the
        bubble-point pressure pb, no gas evolves from the oil and consequently the gas
        solubility remains constant at its maximum value of Rsb.
        - Below the bubble-point pressure, the solution gas is liberated and the
        value of Rs decreases with pressure.

        The following five empirical correlations for estimating the gas solubility are
        given below:

        • Standing’s correlation
        • The Vasquez-Beggs correlation
        • Glaso’s correlation
        • Marhoun’s correlation
        • The Petrosky-Farshad correlation

        """
        sys.path.append(os.path.dirname(__file__))

        method_parts = method.split('_')

        method_class = method_parts[0].capitalize()+''.join(word.capitalize() for word in method_parts[1:])

        # Import the correct method class dynamically
        try:
            module = __import__(method)
            mclass = getattr(module,method_class)
        except (ImportError, AttributeError):
            raise ValueError(f"Method '{method}' not found or invalid.")

        # Create an instance of the class and calculate z-factor
        method_instance = mclass(critical_params,temperature)

        return method_instance(pressures,derivative)

    def rho(pres):
        
        rhor = (rhoSTO+0.01357*Rs*gamma_gas)/fvf

        return rhob*numpy.exp(comp*(pres-pbubble))

    #print(rho(150,2500,65,90,2500,1.23,300,.65,35))

    # def get_density(self):
    # 20. crude oil density
    def rho(T, P, Tsep, Psep, Pb, Bo, Rs, gas_grav, oil_grav):
        """
        The crude oil density is defined as the mass of a unit volume of the
        crude at a specified pressure and temperature. It is usually expressed in
        pounds per cubic foot.

        Function to Calculate Oil Density in lb/ft
        #'T          temperature, °F
        #'P          pressure, psia
        #'Tsep       separator temperature, °F
        #'Psep       separator pressure, psia
        #'Pb         bubble point pressure, psia
        #'Bo         oil formation volume factor, bbl/stb
        #'Rs         solution gas-oil ratio, scf/stb
        #'gas_grav   gas specific gravity
        #'oil_grav   API oil gravity

        """
        oil_grav_sp = 141.5 / (oil_grav + 131.5)
        if (P <= Pb):
            rho_o = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bo)
        else:
            co = oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav)
            Bob = Bo / (math.exp(co * (P - Pb)))
            rho_ob = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bob)
            rho_o = rho_ob * Bo / Bob
        
        return rho_o