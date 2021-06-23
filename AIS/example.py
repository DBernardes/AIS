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
    "image_size": 1024,
}


ais = Artificial_Image_Simulator(
    ccd_operation_mode=dic,
    channel=1,
    gaussian_std=8,
    star_coordinates=[512, 512],
    bias_level=500,
    sparc4_operation_mode="phot",
    image_dir=os.path.join("..", "FITS"),
    star_wavelength_interval=(350, 1150, 50),
    star_temperature=5700,
)

ais.apply_atmosphere_spectral_response()
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response()
ais.create_random_image(1)
