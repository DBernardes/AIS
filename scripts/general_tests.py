# APPLY SERKOWSKI POLARIZATION CURVE

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from AIS.Spectral_Response import Channel


ch = Channel(1)
ch.obj_wavelength = np.asarray([500])
ch.sed = np.zeros((4, 1))
ch.sed[0] = 1
ch.retarder_waveplate_angle = 22.5
ch.retarder_waveplate = "half"

ch._apply_real_waveplate()
print(ch.sed)
