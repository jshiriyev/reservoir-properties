import os, sys

import numpy as np

class Oil():

	# def get_density(self):
    # 20. crude oil density
    def rho(T, P, Tsep, Psep, Pb, Bo, Rs, gas_grav, oil_grav):
        """Function to Calculate Oil Density in lb/ft"""
        #'T          temperature, °F
        #'P          pressure, psia
        #'Tsep       separator temperature, °F
        #'Psep       separator pressure, psia
        #'Pb         bubble point pressure, psia
        #'Bo         oil formation volume factor, bbl/stb
        #'Rs         solution gas-oil ratio, scf/stb
        #'gas_grav   gas specific gravity
        #'oil_grav   API oil gravity
        oil_grav_sp = 141.5 / (oil_grav + 131.5)
        if (P <= Pb):
            rho_o = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bo)
        else:
            co = oil_comp(T, P, Tsep, Psep, Rs, gas_grav, oil_grav)
            Bob = Bo / (math.exp(co * (P - Pb)))
            rho_ob = (350 * oil_grav_sp + 0.0764 * gas_grav * Rs) / (5.615 * Bob)
            rho_o = rho_ob * Bo / Bob
        
        return rho_o

    def rho(pres):
        
        rhor = (rhoSTO+0.01357*Rs*gamma_gas)/fvf

        return rhob*numpy.exp(comp*(pres-pbubble))

    #print(rho(150,2500,65,90,2500,1.23,300,.65,35))

    def fvf(pres,temp,Rs,gamma_gas,gamma_oil):

        CBob = Rs*(gamma_gas/gamma_oil)**0.5+1.25*temp

        Bob = 0.9759+12e-5*CBob**1.2

        return Bob*numpy.exp(co*(pbubble-pres))

    @staticmethod
    def gas_solubility():
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

class Gravity:
    """
    Crude Oil Gravity

    The crude oil density is defined as the mass of a unit volume of the
    crude at a specified pressure and temperature. It is usually expressed in
    pounds per cubic foot. 

    """

    @staticmethod
    def oil(rhoo:float|np.ndarray,rhow:float=62.4):
        """
        The specific gravity of a crude oil is defined as the ratio of the
        density of the oil to that of water. Both densities are measured
        at 60°F and atmospheric pressure:

        Inputs:
        ------
        rhoo = Density of the crude oil, lb/ft3
        rhow = Density of the water, lb/ft3

        Returns:
        -------
        Specific gravity of the oil

        It should be pointed out that the liquid specific gravity is dimensionless,
        but traditionally is given the units 60°/60° to emphasize the fact that
        both densities are measured at standard conditions. The density of the
        water is approximately 62.4 lb/ft3.

        """
        return rhoo/rhow

    @staticmethod
    def API(spgr:float|np.ndarray):
        """
        Calculates API gravity from specific gravity:
        
        Inputs:
        ------
        spgr    : specific gravity of the oil

        Returns:
        -------
        API gravity of the oil

        The API gravities of crude oils usually range from 47° API for the
        lighter crude oils to 10° API for the heavier asphaltic crude oils.

        """
        return 141.5/spgr-131.5

    @staticmethod
    def solg(stock_tank:tuple[float,float],*separators):
        """
        Calculates the specific gravity of the solution gas.

        The specific gravity of the solution gas is described by the weighted
        average of the specific gravities of the separated gas from each separator.

        This weighted-average approach is based on the separator gas-oil ratio.
        
        The stock_tank is a tuple with Rst and gst values:
            Rst = gas-oil ratio from the stock tank, scf/ STB
            gst = gas gravity from the stock tank
        
        For each separator Rsep and gsep values must be provided:
            Rsep = separator gas-oil ratio, scf/STB
            gsep = separator gas gravity
        
        Returns:
        -------
        Specific gravity of the solution gas.

        """
        Rst,gst = stock_tank

        upper = Rst*gst
        lower = Rst

        for sep in separators:
            Rsep,gsep = sep

            upper += Rsep*gsep
            lower += Rsep

        return upper/lower

    @staticmethod
    def gas_correction(Tsep,psep,API,gg):
        """Method to calculate corrected gas gravity:

        Realizing that the value of the specific gravity of the gas depends on
        the conditions under which it is separated from the oil, Vasquez and
        Beggs proposed that the value of the gas specific gravity as obtained
        from a separator pressure of 100 psig be used in the above equation. This
        reference pressure was chosen because it represents the average field
        separator conditions.

        Tsep    actual separator temperature, °F
        psep    actual separator pressure, psia
        API     API oil gravity
        gg      gas gravity at the actual separator conditions of psep and Tsep

        The gas gravity used to develop all the correlations reported by the
        authors was that which would result from a two-stage separation. The
        first-stage pressure was chosen as 100 psig and the second stage was the
        stock tank. If the separator conditions are unknown, the unadjusted gas
        gravity may be used.

        Returns:
        -------
        Gas gravity at the reference separator pressure

        """
        return gg*(1+5.912e-5*API*Tsep*np.log10(psep/(100.+14.7)))

if __name__ == "__main__":

    sg = spgr(53)

    print(sg)
    print(API(sg))

    print(solg((58,1.296),(724,0.743),(202,0.956)))