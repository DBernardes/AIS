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
wv = np.linspace(400, 1100, 10)
sed = np.zeros((4, 10))
sed[0, :] = 1

ch.sed = sed
ch.obj_wavelength = wv
ch.write_sparc4_operation_mode("polarimetry", retarder_waveplate_angle=22.5)
ch._apply_real_polarizer(np.ones(10))

print(ch.sed)
