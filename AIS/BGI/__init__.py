"""Back Ground Image Generator.


This package calculates a back ground image that simulates and bias image of
the SPARC4 cameras
"""


from .BGI import Back_Ground_Image
from FC import Flux_Calculation
from Telescope_SR import Telescope_Spectral_Response
from Atmosphere_SR import Atmosphere_Spectral_Response
from S4_SR import (Abstract_SPARC4_Spectral_Response,
                   Concrete_SPARC4_Spectral_Response_1,
                   Concrete_SPARC4_Spectral_Response_2,
                   Concrete_SPARC4_Spectral_Response_3,
                   Concrete_SPARC4_Spectral_Response_4)
