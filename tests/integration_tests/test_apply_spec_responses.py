# # -*- coding: utf-8 -*-


import unittest
from math import e

import astropy.io.fits as fits
import numpy as np
import pytest
from numpy import ndarray
from scipy.constants import c, h, k
from scipy.interpolate import (
    PchipInterpolator,
    UnivariateSpline,
    interp1d,
    splev,
    splrep,
)

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator


class Test_Apply_Spectral_Response(unittest.TestCase):
    channel = 1
    ccd_temp = -70
    ccd_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "readout": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 200,
    }
    mag = 12
    wv_range = (350, 1100, 10)
    star_temp = 5700

    @classmethod
    def setUp(cls):
        # When using SetUp() rather than SetUpClass(), this method will be
        # executed for each new test
        cls.sed = np.zeros((4, 10))
        cls.sed[0] = np.ones(10)
        cls.wv = np.linspace(*cls.wv_range)
        cls.ais = Artificial_Image_Simulator(cls.ccd_mode, cls.channel, cls.ccd_temp)
        cls.ais.write_source_sed(cls.wv, cls.sed)
        cls.ais.create_sky_sed("new")

    def _interpolate(self, wavelength, spectral_response) -> ndarray:
        b = PchipInterpolator(wavelength, spectral_response)
        spectral_response = b(self.wv)

        return spectral_response

    # -----------------------------------------------------------------------------

    def test_write_source_sed(self):
        assert np.allclose(self.AIS.source_sed, self.sed)
        assert np.allclose(self.AIS.wavelength, self.wv)

    def test_apply_telescope(self):
        self.AIS.apply_telescope_spectral_response()

    def test_interpolate(self):
        spec_response = self._interpolate(self.wv, self.sed[0])
        assert np.mean(spec_response) == 1
