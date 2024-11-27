# APPLY SERKOWSKI POLARIZATION CURVE

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

sed = np.linspace(1, 100, 100)
wavelength = np.linspace(350, 1100, 100)
ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}
ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
ais.write_source_sed(wavelength, sed)
ais.apply_linear_polarization(50, 0)
