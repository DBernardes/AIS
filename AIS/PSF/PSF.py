"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class calculates the star flux distribution
based on a gaussian 2D distributiopn
"""

from FC import Flux_Calculation
from Telescope_SR import Telescope_Spectral_Response
from Atmosphere_SR import Atmosphere_Spectral_Response
from CHC import Abstract_Channel_Creator


class Point_Spread_Function:
    """Point Spread Function Class."""

    def __init__(self, Abstract_Channel_Creator, ccd_gain):
        self.CHC = Abstract_Channel_Creator
        self.FC = Flux_Calculation
        self.TSR = Telescope_Spectral_Response
        self.ASR = Atmosphere_Spectral_Response
        self.ccd_gain = ccd_gain
