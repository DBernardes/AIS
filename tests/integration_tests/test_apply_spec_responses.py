# # -*- coding: utf-8 -*-


import unittest
from math import e

import astropy.io.fits as fits
import numpy as np
import pytest
from scipy.constants import c, h, k

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
        cls.SED = np.zeros((4, 10))
        cls.SED[0] = np.ones(10)
        cls.WV = np.linspace(*cls.wv_range)
        cls.AIS = Artificial_Image_Simulator(cls.ccd_mode, cls.channel, cls.ccd_temp)
        # cls.AIS.write_source_sed(cls.WV, cls.SED)
        cls.AIS.create_source_sed("blackbody", cls.mag, cls.wv_range, cls.star_temp)
        cls.AIS.create_sky_sed("new")

    def black_body(wavelength, temperature):
        numerator = 2 * h * c**2
        denominator = (wavelength**5) * (
            np.exp((h * c) / (wavelength * k * temperature)) - 1
        )

        return numerator / denominator

    # -----------------------------------------------------------------------------

    def test_apply_telescope(self):
        self.AIS.apply_telescope_spectral_response()
        print("\n", self.AIS.source_sed[0])

    def test_create_blackbody(self):
        pass
