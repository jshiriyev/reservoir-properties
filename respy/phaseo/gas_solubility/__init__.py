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
crude oil, is shown in Figure 2-7. As the pressure is reduced
from the initial reservoir pressure pi, to the bubble-point pressure pb, no
gas evolves from the oil and consequently the gas solubility remains constant
at its maximum value of Rsb. Below the bubble-point pressure, the
solution gas is liberated and the value of Rs decreases with pressure. The
following five empirical correlations for estimating the gas solubility are
given below:

• Standing’s correlation
• The Vasquez-Beggs correlation
• Glaso’s correlation
• Marhoun’s correlation
• The Petrosky-Farshad correlation

"""
# from ._standings_correlation import Standing
# from ._vasquez_beggs_correlation import VasquezBeggs
# from ._glasos_correlation import Glaso
# from ._marhouns_correlation import Marhoun
# from ._petrosky_farshad_correlation import PetroskyFarshad