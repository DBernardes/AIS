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

from sbpy.calib import Vega, vega_fluxd

print(vega_fluxd.get()["Johnson V"].value)
