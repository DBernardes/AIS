# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import os
from math import cos

import numpy as np

from AIS.Spectral_Response import Channel, _utils

ch = Channel(1)
wv = np.linspace(400, 1100, 100)
sed = np.zeros((4, 100))
sed[0, :] = np.ones(100) * 5
sed[1, :] = np.ones(100) * 5

ch.sed = sed
ch.obj_wavelength = wv
ch.write_sparc4_operation_mode("polarimetry", retarder_waveplate_angle=180)
# ch._apply_real_waveplate()
# # ch._apply_analyzer()
# print(ch.sed)
