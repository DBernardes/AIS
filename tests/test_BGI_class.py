# -*- coding: utf-8 -*-
"""Tests of the Backgorund Image class.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""


import unittest

import numpy as np
from photutils.datasets import make_noise_image

from AIS.Background_Image import Background_Image
from AIS.Noise import Noise


class Test_Back_Ground_Image(unittest.TestCase):
    CCD_OP_MODE = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "readout": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 1024,
    }
    TEMP = -70
    BGI_CONV = Background_Image(CCD_OP_MODE, 1, TEMP)
    CCD_GAIN = BGI_CONV.ccd_gain
    PIXEL_SENSIBILITY = 0.03
    SEED = 5

    @classmethod
    def setUpClass(cls):
        noise = Noise(1)
        cls.read_noise = noise.calculate_read_noise(cls.CCD_OP_MODE)
        cls.dark_noise = (
            noise.calculate_dark_current(cls.TEMP) * cls.CCD_OP_MODE["t_exp"]
        )

    # ------------------------ Initialize the class --------------------------

    def test_ccd_operation_mode(self):
        assert self.BGI_CONV.ccd_operation_mode == self.CCD_OP_MODE

    def test_bias_level(self):
        assert self.BGI_CONV.bias_level == 500

    def test_channel(self):
        assert self.BGI_CONV.channel == 1

    def test_noise_factor(self):
        assert self.BGI_CONV._NOISE_FACTOR == 1

    def test_read_noise(self):
        assert self.BGI_CONV.read_noise == 6.67

    def test_dark_noise(self):
        assert (
            self.BGI_CONV.dark_noise
            == 5.8597559895090484e-05 * self.CCD_OP_MODE["t_exp"]
        )

    def test_get_ccd_gain(self):
        assert self.BGI_CONV.ccd_gain == 3.37

    def test_create_bias_background(self):
        bbg = self.BGI_CONV.create_bias_background(self.SEED)

        image_size = self.CCD_OP_MODE["image_size"]
        shape = (image_size, image_size)
        noise_adu = self.read_noise / self.CCD_GAIN
        bias_background = make_noise_image(
            shape, distribution="gaussian", mean=500, stddev=noise_adu, seed=self.SEED
        )

        assert np.allclose(bbg, bias_background)

    def test_create_dark_background(self):
        image_size = self.CCD_OP_MODE["image_size"]
        em_gain = self.CCD_OP_MODE["em_gain"]
        binn = self.CCD_OP_MODE["binn"]
        shape = (image_size, image_size)
        dark_level = 500 + self.dark_noise * em_gain * binn**2 / self.CCD_GAIN

        noise = (
            np.sqrt(self.read_noise**2 + self.dark_noise * em_gain**2 * binn**2)
            / self.CCD_GAIN
        )

        dark_background = make_noise_image(
            shape,
            distribution="gaussian",
            mean=dark_level,
            stddev=noise,
            seed=self.SEED,
        )

        assert np.allclose(
            dark_background, self.BGI_CONV.create_dark_background(self.SEED)
        )

    def test_create_flat_image(self):
        em_gain = self.CCD_OP_MODE["em_gain"]
        binn = self.CCD_OP_MODE["binn"]
        image_size = self.CCD_OP_MODE["image_size"]
        FLAT_LEVEL = 2**14

        poisson_noise = FLAT_LEVEL / self.CCD_GAIN
        shape = (image_size, image_size)

        noise = (
            np.sqrt(
                self.read_noise**2
                + (poisson_noise * (1 + self.PIXEL_SENSIBILITY) + self.dark_noise)
                * em_gain**2
                * binn**2
            )
            / self.CCD_GAIN
        )

        flat_background = make_noise_image(
            shape,
            distribution="gaussian",
            mean=FLAT_LEVEL,
            stddev=noise,
            seed=self.SEED,
        )

        assert np.allclose(
            flat_background, self.BGI_CONV.create_flat_background(self.SEED)
        )

    def test_create_sky_background(self):
        t_exp = self.CCD_OP_MODE["t_exp"]
        em_gain = self.CCD_OP_MODE["em_gain"]
        binn = self.CCD_OP_MODE["binn"]
        image_size = self.CCD_OP_MODE["image_size"]

        sky_background = (
            500 + (self.dark_noise + 10) * t_exp * em_gain * binn**2 / self.CCD_GAIN
        )

        noise = (
            np.sqrt(
                self.read_noise**2
                + (10 * t_exp + self.dark_noise) * em_gain**2 * binn**2
            )
            / self.CCD_GAIN
        )

        shape = (image_size, image_size)
        sky_background = make_noise_image(
            shape,
            distribution="gaussian",
            mean=sky_background,
            stddev=noise,
            seed=self.SEED,
        )

        assert np.allclose(
            sky_background, self.BGI_CONV.create_sky_background(10, self.SEED)
        )
