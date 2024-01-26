import os
import sys

import numpy as np
import pandas as pd

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}
dict = {"angle": [], "ord": [], "extra": []}
channel = 1
ais = Artificial_Image_Simulator(
    ccd_operation_mode, channel_id=channel, ccd_temperature=-70
)
for i in range(1):
    ais.create_source_sed(
        calculation_method="blackbody",
        magnitude=15,
        wavelength_interval=(400, 1100, 100),
        temperature=5700,
    )
    ais.apply_linear_polarization(100, 0)
    ais.create_sky_sed(moon_phase="new")
    ais.apply_sparc4_spectral_response(
        "polarimetry",
        retarder_waveplate="half",
        retarder_waveplate_angle=i * 22.5,
    )
    ais._integrate_sed()
