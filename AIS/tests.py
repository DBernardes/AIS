# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import matplotlib.pyplot as plt
import numpy as np

from Channel_Creator import Abstract_Channel_Creator

l_init, l_final, l_step = 350, 1100, 50
wl = range(l_init, l_final, l_step)
n = len(wl)
specific_flux = np.ones((4, n))
abs = Abstract_Channel_Creator(-70, "pol")
abs.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
