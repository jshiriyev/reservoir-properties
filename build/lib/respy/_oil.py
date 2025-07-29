import importlib
import os
import sys

import numpy as np

from ._fluid import Fluid

from .phaseo._crude_oil_system import CrudeOilSystem

class OilPhase(CrudeOilSystem):

    def __init__(self,sgsg:float,gAPI:float,temp:float):
        """
        Initialize a phaseo instance.

        Args:
        ----
        sgsg (float): Specific gravity of solution gas.
        gAPI (float): API gravity of oil.
        temp (float): Temperature in Fahrenheit.

        """
        self._sgsg  = sgsg
        self._gAPI  = gAPI
        self._temp  = temp

    @property
    def sgsg(self):
        """Get the specific gravity of the solution gas."""
        return self._sgsg

    @property
    def gAPI(self):
        """Get the API gravity of the oil."""
        return self._gAPI

    @property
    def temp(self):
        """Get the temperature of the reservoir."""
        return self._temp

    def __call__(self,bpp:float=None,Rsb:float=None,method="standings_correlation",**kwargs):

        sys.path.append(os.path.dirname(__file__))

        method_parts = method.split('_')
        method_class = method_parts[0].capitalize()+''.join(word.capitalize() for word in method_parts[1:])
        
        # Import the correct method class dynamically
        try:
            module = importlib.import_module(f"phaseo._{method}")
            self.__corr = getattr(module,method_class)
        except (ImportError, AttributeError):
            raise ValueError(f"Method '{method}' not found or invalid.")

        if bpp is None:
            self._bpp = self.__corr.bpp(Rsb,self.sgsg,self.gAPI,self.temp,**kwargs)
        else:
            self._bpp = bpp

        return self

    @property
    def bpp(self):
        """Get the bubble point pressure."""
        return self._bpp
    
    def gass(self,p:np.ndarray,**kwargs):
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

        • Standing’s correlation (standing)
        • The Vasquez-Beggs correlation (vasquez_beggs)
        • Glaso’s correlation (glaso)
        • Marhoun’s correlation (marhoun)
        • The Petrosky-Farshad correlation (petrosky_farshad)

        """
        return self.__corr.gass(p,self.bpp,self.sgsg,self.gAPI,self.temp,**kwargs)

    def fvf(self,p:np.ndarray,**kwargs):
        """
        The oil formation volume factor, Bo, is defined as the ratio of the volume
        of oil (plus the gas in solution) at the prevailing reservoir temperature
        and pressure to the volume of oil at standard conditions. Bo is always
        greater than or equal to unity. The oil formation volume factor can be
        expressed mathematically as:

        Bo = (Vo)p,T / (Vo)sc

        where:

        Bo      = oil formation volume factor, bbl/STB
        (Vo)p,T = volume of oil under reservoir pressure p and temperature T, bbl
        (Vo)sc  = volume of oil is measured under standard conditions, STB

        A typical oil formation factor curve, as a function of pressure for an
        undersaturated crude oil (pi > pb), is as following:

         - As the pressure is reduced below the initial reservoir pressure pi,
        the oil volume increases due to the oil expansion. This behavior results
        in an increase in the oil formation volume factor and will continue until
        the bubble-point pressure is reached.

         - At pb, the oil reaches its maximum expansion and consequently
        attains a maximum value of Bob for the oil formation volume factor.

        - As the pressure is reduced below pb, volume of the oil and Bo are
        decreased as the solution gas is liberated. When the pressure is reduced
        to atmospheric pressure and the temperature to 60°F, the value of Bo is
        equal to one.

        Most of the published empirical Bo correlations utilize the following
        generalized relationship:

        Bo = f (Rs,gg,go,T)

        Six different methods of predicting the oil formation volume factor are
        presented below:

        • Standing’s correlation
        • The Vasquez-Beggs correlation
        • Glaso’s correlation
        • Marhoun’s correlation
        • The Petrosky-Farshad correlation
        • Other correlations

        It should be noted that all the correlations could be used for any pressure
        equal to or below the bubble-point pressure.

        """
        return self.__corr.fvf(p,self.bpp,self.sgsg,self.gAPI,self.temp,**kwargs)

    def rho(self,p,Bo,Rs,sgsg,gAPI,temp,Tsep,psep):
        """
        The crude oil density is defined as the mass of a unit volume of the
        crude at a specified pressure and temperature. It is usually expressed in
        pounds per cubic foot.

        Function to Calculate Oil Density in lb/ft3

        T     : temperature, °F
        P     : pressure, psia
        sgsg  : gas specific gravity
        gAPI  : API oil gravity

        Tsep  : separator temperature, °F
        psep  : separator pressure, psia
        Pb    : bubble point pressure, psia
        Bo    : oil formation volume factor, bbl/stb
        Rs    : solution gas-oil ratio, scf/stb

        """
        # rhor = (rhoSTO+0.01357*Rs*gamma_gas)/fvf
        # return rhob*numpy.exp(comp*(pres-pbubble))
        #print(rho(150,2500,65,90,2500,1.23,300,.65,35))

        API_sp = 141.5 / (gAPI + 131.5)

        if (P <= Pb):
            rho_o = (350 * API_sp + 0.0764 * sgsg * Rs) / (5.615 * Bo)
        else:
            co = oil_comp(T, P, Tsep, psep, Rs, sgsg, gAPI)
            Bob = Bo / (math.exp(co * (P - Pb)))
            rho_ob = (350 * API_sp + 0.0764 * sgsg * Rs) / (5.615 * Bob)
            rho_o = rho_ob * Bo / Bob
        
        return rho_o

    def comp(self,p:np.ndarray,**kwargs):
        """
        Isothermal compressibility coefficients are required in solving many
        reservoir engineering problems, including transient fluid flow problems,
        and they are also required in the determination of the physical properties
        of the undersaturated crude oil.

        By definition, the isothermal compressibility of a substance is defined
        mathematically by the following expression:

        c = -1 / V * (\\partial V / \\partial p)_T

        For a crude oil system, the isothermal compressibility coefficient of
        the oil phase co is defined for pressures above the bubble-point by one of
        the following equivalent expressions:

        There are several correlations that are developed to estimate the oil compressibility
        at pressures above the bubble-point pressure, i.e., undersaturated
        crude oil system. Three of these correlations are presented below:

        • The Vasquez-Beggs correlation
        • The Petrosky-Farshad correlation
        • McCain’s correlation

        """
        return self.__corr.comp(p,*args,**kwargs)

    def visc(self):

        pass