# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import os

from Artificial_Image_Simulator import Artificial_Image_Simulator

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 200,
}


ais = Artificial_Image_Simulator(
    ccd_operation_mode=dic,
    channel=1,
    star_coordinates=[100, 100],
    bias_level=500,
    sparc4_operation_mode="pol",
    image_dir=os.path.join("..", "FITS"),
    wavelength_interval=(400, 1150, 50),
    star_temperature=5700,
    star_magnitude=18,
)

# ais.apply_atmosphere_spectral_response()
# ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response()
# ais.create_artificial_image()
# ais.create_background_image()
# ais.create_bias_image()
# ais.create_dark_image()
# ais.create_flat_image()
ais.create_random_image(10)
