# # -*- coding: utf-8 -*-
# Tests of the Sky Class
# Oct 24th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest, unittest
from AIS.Spectral_Energy_Distribution import Sky
from scipy.interpolate import splev, splrep, interp1d
from scipy.constants import c, h
from sbpy.calib import vega_fluxd


class Test_Sky(unittest.TestCase):
    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    BASE_PATH = os.path.join("AIS", "Spectral_Energy_Distribution")
    CSV_FILE = "moon_magnitude.csv"

    @classmethod
    def setUpClass(cls):
        cls.sky = Sky()

    def test_read_csv(self):
        file_name = os.path.join(self.BASE_PATH, self.CSV_FILE)
        ss = pd.read_csv(file_name)
        wavelenght = ss["wavelength"]
        value = ss["full"]

        tup = self.sky._read_csv(self.CSV_FILE, "full")
        assert np.allclose(tup[0], wavelenght)
        assert np.allclose(tup[1], value)

    def test_calculate_sed(self):
        wavelength, mags = self.sky._read_csv(self.CSV_FILE, "full")
        sed = self.sky._calculate_photons_density(mags)
        temp = self.sky._interpolate_spectral_distribution(
            wavelength, sed, self.OBJ_WAVELENGTH
        )
        sed = np.zeros((4, self.OBJ_WAVELENGTH.shape[0]))
        sed[0] = temp

        new_sed = self.sky.calculate_sed("full", self.OBJ_WAVELENGTH)
        assert np.allclose(sed, new_sed)
