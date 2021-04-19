"""
Background Image Class.
=======================

This is the Background Image Class used to generate a back ground image like
a bias image produced by the SPARC4 cameras.
"""

from FC import Flux_Calculation
from Telescope_SR import Telescope_Spectral_Response
from Atmosphere_SR import Atmosphere_Spectral_Response
from CHC import Abstract_Channel_Creator


class Background_Image:
    """Background Image Class."""

    def __init__(self, Abstract_Channel_Creator, ccd_gain):
        self.CHC = Abstract_Channel_Creator
        self.FC = Flux_Calculation
        self.TSR = Telescope_Spectral_Response
        self.ASR = Atmosphere_Spectral_Response
        self.ccd_gain = ccd_gain
