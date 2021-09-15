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
from Telescope_Spectral_Response import Telescope_Spectral_Response

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
specific_flux = np.ones((4, n))

tel_sr = Telescope_Spectral_Response()
flux = tel_sr.apply_telescope_spectral_response(specific_flux, l_init, l_final, l_step)
