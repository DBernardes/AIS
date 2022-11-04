# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

This is a script for the development of general tests

@author: denis
"""

# -----------------------------------------------------------------------------------------
# Test the PSF in hte polarimetric mode
# from AIS.Point_Spread_Function import Point_Spread_Function
# from tests.AIS_spectral_response_curves import ccd_operation_mode
# import matplotlib.pyplot as plt
# psf = Point_Spread_Function(ccd_operation_mode, 1)
# image = psf.create_star_image((50, 50), 100, 100)
# plt.imshow(image)
# plt.show()
# -----------------------------------------------------------------------------------------
# from AIS.Header import Header
# Pint the header
# from tests.AIS_spectral_response_curves import ccd_operation_mode
# hdr = Header(ccd_operation_mode)
# header = hdr.create_header()
# print(repr(header))
# -----------------------------------------------------------------------------------------
# Test the blackbody profile
# from AIS.Spectral_Energy_Distribution import Source
# import matplotlib.pyplot as plt
# import numpy as np
# source = Source()
# sed = source.calculate_sed(
#     'blackbody', (350, 1100, 100), 1, 1, temperature=5700)
# plt.plot(np.linspace(350, 1100, 100), sed)
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the sbpy package
# from sbpy.calib import Vega, vega_fluxd

# print(vega_fluxd.get()["Johnson V"].value)
# -----------------------------------------------------------------------------------------
# Test the _read_spectral_library function
# from AIS.Spectral_Energy_Distribution import Source
# from matplotlib import pyplot as plt
# source = Source()
# wv, sed = source.calculate_sed(
#     'spectral_library', 10, spectral_type='G')
# plt.plot(wv, sed)
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the Sky Class
# from AIS.Spectral_Energy_Distribution import Sky
# import matplotlib.pyplot as plt
# from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength
# sky = Sky()
# sed = sky.calculate_sed('new', obj_wavelength)
# plt.plot(obj_wavelength, sed)
# plt.show()
# -----------------------------------------------------------------------------------------
from AIS.Spectral_Energy_Distribution import Source
from matplotlib import pyplot as plt

source = Source()
sv, sed = source.calculate_sed('spectral_library', 10, spectral_type='aaa')
plt.plot(sv, sed)
plt.show()
