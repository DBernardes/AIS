# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import os
from copy import copy
from math import cos, sin

import matplotlib.pyplot as plt
import numpy as np

from AIS.Spectral_Response import Channel, _utils

ch = Channel(1)
wv = np.linspace(400, 1100, 100)
sed = np.zeros((4, 100))
sed[0, :] = np.ones(100)
sed[1, :] = np.ones(100) * 5
sed[2, :] = np.ones(100) * 5

ch.sed = sed
ch.obj_wavelength = wv
ch.write_sparc4_operation_mode("polarimetry", retarder_waveplate_angle=22.5)
ch._apply_ideal_waveplate()
ch._apply_analyzer()
plt.plot(ch.sed[0, :])
plt.plot(ch.sed[1, :])
plt.show()
