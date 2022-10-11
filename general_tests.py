# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

This is a script for the development of general tests

@author: denis
"""

# -----------------------------------------------------------------------------------------
from AIS.Point_Spread_Function import Point_Spread_Function
from tests.AIS_spectral_response_curves import ccd_operation_mode
import matplotlib.pyplot as plt
psf = Point_Spread_Function(ccd_operation_mode, 1)
image = psf.create_star_image((50, 50), 100, 100)
plt.imshow(image)
plt.show()
# -----------------------------------------------------------------------------------------
