from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
import pandas as pd
import os
import numpy as np
from scipy.constants import c, h
import matplotlib.pyplot as plt
from sys import exit

ccd_operation_mode = {
    'em_mode': 'Conv',
    'em_gain': 1,
    'preamp': 1,
    'readout': 1,
    'binn': 1,
    't_exp': 1,
    'image_size': 100
}

ais = Artificial_Image_Simulator(
    ccd_operation_mode, channel_id=4, ccd_temperature=-70)
ais.create_source_sed('blackbody', 12, (400, 1100, 100), 5700)
ais.create_sky_sed(moon_phase='first quarter')
ais.apply_atmosphere_spectral_response(air_mass=1.019)
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response(
    acquisition_mode='photometry')
ais.create_background_image('tests/fits')
