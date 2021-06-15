# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Channel_Creator import Abstract_Channel_Creator

wavelength_interv = range(350, 1100, 50)
n = len(wavelength_interv)
specific_flux = np.ones((4, n))
chc = Abstract_Channel_Creator(sparc4_operation_mode="pol")
ord, extra_ord = chc.apply_sparc4_spectral_response(specific_flux, 350, 1100, 50)
print(ord, extra_ord)
