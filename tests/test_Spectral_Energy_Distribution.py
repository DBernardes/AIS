# # -*- coding: utf-8 -*-
# Tests of the Spectral Energy Distributiion Class
# Oct 24th 2022
# @author: denis
#

import os
import unittest

import numpy as np
import pandas as pd
from sbpy.calib import vega_fluxd
from scipy.constants import c, h, k
from scipy.interpolate import splev, splrep

from AIS.Spectral_Energy_Distribution import Spectral_Energy_Distribution


class Test_Spectral_Energy_Distribution(unittest.TestCase):
    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.linspace(100, 1000, 100)
    WV = np.linspace(350, 1100, 100)
    SPEC_RESPONSE = np.ones(100)
    BASE_PATH = os.path.join("AIS", "Spectral_Energy_Distribution")

    @classmethod
    def setUpClass(cls):
        cls.sed = Spectral_Energy_Distribution()
        spl = splrep(cls.WV, cls.SPEC_RESPONSE)
        cls.interpolated_sed = splev(cls.WV, spl)

    def test_interpolate(self):
        interpolated_sed = self.sed._interpolate_spectral_distribution(
            self.WV, self.SPEC_RESPONSE, self.OBJ_WAVELENGTH
        )
        assert np.allclose(self.interpolated_sed, interpolated_sed)

    def test_calc_sed(self):
        assert np.allclose(self.sed.calculate_sed(), self.SED)

    def test_calculate_photons_density(self):
        magnitude = 10
        TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
        EFFECT_WAVELENGTH = 545  # nm
        S_0 = 3.631e-2  # W/m2/m
        photons_density = (
            S_0
            * 10 ** (-magnitude / 2.5)
            * TELESCOPE_EFFECTIVE_AREA
            * EFFECT_WAVELENGTH
            * 1e-9
            / (h * c)
        )
        assert self.sed._calculate_photons_density(magnitude) == photons_density
