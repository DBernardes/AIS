# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import os
from random import gauss

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Artificial_Image_Simulator import Artificial_Image_Simulator
from SPARC4_Spectral_Response import Concrete_SPARC4_Spectral_Response_1
from Spectrum_Calculation import Spectrum_Calculation
from Telescope_Spectral_Response import Telescope_Spectral_Response
from tests.SPARC4_SR_curves import *
from tests.test_AIS_operation import multiply_matrices
from tests.test_S4_SR_class import c1_s4_sr

# dic = {
#     "em_mode": 0,
#     "em_gain": 1,
#     "preamp": 1,
#     "hss": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "ccd_temp": -70,
#     "image_size": 200,
# }

# gaussian_std = 8
# ais = Artificial_Image_Simulator(
#     dic,
#     image_dir=os.path.join("..", "FITS"),
#     sparc4_operation_mode="pol",
#     gaussian_std=gaussian_std,
# )
# ais.apply_sparc4_spectral_response()

# -------------------------------------------------------------------------------

l_init, l_final, l_step = 350, 1150, 50
wavelength_interval = range(l_init, l_final, l_step)
n = len(wavelength_interval)
sc = Spectrum_Calculation(5700, l_init, l_final, l_step)
specific_flux = sc.calculate_sky_specific_flux()

c1_s4sr = Concrete_SPARC4_Spectral_Response_1()
c1_s4sr.write_specific_flux(specific_flux.copy(), wavelength_interval)

c1_s4sr.apply_calibration_wheel()
specific_flux = multiply_matrices(calibration_wheel, specific_flux)

c1_s4sr.apply_retarder()
specific_flux = multiply_matrices(retarder, specific_flux)

c1_s4sr.apply_analyser()
specific_flux = multiply_matrices(analyser_ordinary_ray, specific_flux)

c1_s4sr.apply_collimator()
specific_flux = np.multiply(colimator_transmitance, specific_flux[0, :])

c1_s4sr.apply_dichroic()
specific_flux = np.multiply(dichroic_c1, specific_flux)

c1_s4sr.apply_camera()
specific_flux - np.multiply(camera_c1, specific_flux)

c1_s4sr.apply_ccd()
specific_flux = np.multiply(ccd_transmitance_c1, specific_flux)
print(np.allclose(specific_flux, c1_s4sr.specific_ordinary_ray))
