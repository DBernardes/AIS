import os

import matplotlib.pyplot as plt
import numpy as np

from AIS.Spectral_Response import Channel
from AIS.Spectral_Response._utils import calculate_retarder_matrix

wv = np.linspace(400, 1100, 20)
sed = np.zeros((4, 20))
sed[0, :] = 1
sed[1, :] = 1

ch = Channel(1)
ch.write_sparc4_operation_mode(
    "polarimetry", retarder_waveplate="quarter", retarder_waveplate_angle=90
)
ch.sed = sed
ch.obj_wavelength = wv
ch._apply_real_waveplate()
print(ch.sed)
