# -*- coding: utf-8 -*-

"""
Flux Calculation class tests.

This script tests the operation of the Flux Calculation Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import os
import unittest
from math import pi

import numpy as np
import pytest
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image, make_noise_image

from AIS.Point_Spread_Function import Point_Spread_Function


class Test_PSF(unittest.TestCase):
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
    PSF = Point_Spread_Function(CCD_OP_MODE, 1)
    SPARC4_POL_SEPARATION = 20  # pix
    SPARC4_PLATE_SCALE = 0.35  # arcsec/pix
    SPREADSHEET_PATH = os.path.join("AIS", "Point_Spread_Function", "preamp_gains.csv")
    SEEING = 1
    SEED = 5

    @classmethod
    def setUpClass(cls):
        pass

    # ------------------------ Initialize the class --------------------------

    def test_ccd_operation_mode(self):
        assert self.PSF.ccd_operation_mode == self.CCD_OP_MODE

    def test_channel(self):
        assert self.PSF.channel == 1

    def test_get_ccd_gain(self):
        assert self.PSF.ccd_gain == 3.37

    # ------------------------------------------------------------------------------------------------------

    def test_create_table(self):
        image_size = self.CCD_OP_MODE["image_size"]
        half_img_size = image_size // 2
        star_coordinates = (half_img_size, half_img_size)
        gaussian_std = self.SEEING / self.SPARC4_PLATE_SCALE
        binn = self.CCD_OP_MODE["binn"]
        x_coord = half_img_size
        y_coord = half_img_size
        table = Table()
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std / binn / 2]
        table["y_stddev"] = [gaussian_std / binn / 2]
        table["theta"] = np.radians(np.array([0]))
        self.PSF._create_table(star_coordinates, self.SEEING)
        assert self.PSF.table == table

    def test_calculate_npix_star(self):
        gaussian_std = self.SEEING / (
            self.SPARC4_PLATE_SCALE * self.CCD_OP_MODE["binn"] * 2
        )
        fwhm = 2.355 * gaussian_std
        psf_star = 3 * fwhm
        npix = pi * psf_star**2
        assert self.PSF.calculate_npix_star(self.SEEING) == npix

    def test_calculate_gaussian_amplitude(self):
        photons_per_second = 100
        t_exp = self.CCD_OP_MODE["t_exp"]
        em_gain = self.CCD_OP_MODE["em_gain"]
        binn = self.CCD_OP_MODE["binn"]
        self.PSF._create_table((100, 100), self.SEEING)

        gaussian_amplitude = (
            photons_per_second
            * t_exp
            * em_gain
            * binn**2
            / (
                self.PSF.ccd_gain
                * 2
                * pi
                * self.PSF.table["x_stddev"]
                * self.PSF.table["y_stddev"]
            )
        )

        assert gaussian_amplitude == self.PSF._calculate_gaussian_amplitude(
            photons_per_second
        )

    def test_make_noise_image(self):
        self.PSF._create_table((50, 50), self.SEEING)
        gaussian_amplitude = self.PSF._calculate_gaussian_amplitude(100)
        self.PSF.table["amplitude"] = gaussian_amplitude
        self.PSF.seed = self.SEED
        image_noise = self.PSF._make_noise_image()

        table = self.PSF.table
        new_image = (
            make_noise_image(
                self.PSF.shape, "poisson", gaussian_amplitude, seed=self.SEED
            )
            - gaussian_amplitude
        )
        table["amplitude"] = [1]
        new_image *= make_gaussian_sources_image(self.PSF.shape, table)
        assert np.allclose(new_image, image_noise)

    def test_create_image_ordinary_ray(self):
        self.PSF._create_table((100, 100), self.SEEING)
        gaussian_amplitude = self.PSF._calculate_gaussian_amplitude(100)
        self.PSF.table["amplitude"] = gaussian_amplitude
        self.PSF.seed = self.SEED
        star_image = make_gaussian_sources_image(self.PSF.shape, self.PSF.table)
        star_image += self.PSF._make_noise_image()

        new_image = self.PSF._create_image_ordinary_ray(100)

        assert np.allclose(star_image, new_image)

    def test_create_image_extra_ordinary_ray(self):
        self.PSF._create_table((100, 100), self.SEEING)
        self.PSF.seed = self.SEED
        new_image = self.PSF._create_image_extra_ordinary_ray(100)

        self.PSF._create_table((100, 100), self.SEEING)
        gaussian_amplitude = self.PSF._calculate_gaussian_amplitude(100)
        self.PSF.table["amplitude"] = gaussian_amplitude
        self.PSF.table["x_mean"] -= self.SPARC4_POL_SEPARATION
        self.PSF.table["y_mean"] -= self.SPARC4_POL_SEPARATION
        star_image = make_gaussian_sources_image(self.PSF.shape, self.PSF.table)
        star_image += self.PSF._make_noise_image()

        assert np.allclose(star_image, new_image)

    def test_creat_star_image(self):
        self.PSF.seed = self.SEED
        self.PSF._create_table((100, 100), self.SEEING)
        gaussian_amplitude = self.PSF._calculate_gaussian_amplitude(100)
        star_image = self.PSF._create_image_ordinary_ray(100)
        star_image += self.PSF._create_image_extra_ordinary_ray(100)

        new_image = self.PSF.create_star_image(
            (100, 100), 100, 100, self.SEEING, self.SEED
        )

        assert np.allclose(star_image, new_image, atol=5 * np.sqrt(gaussian_amplitude))
