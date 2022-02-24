# -*- coding: utf-8 -*-

"""
Flux Calculation class tests.

This script tests the operation of the Flux Calculation Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import matplotlib.pyplot as plt
import numpy as np
import pytest
from AIS.Channel_Creator import Concrete_Channel_1
from AIS.Point_Spread_Function import Point_Spread_Function
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image, make_noise_image

from .AIS_spectral_response_curves import ccd_operation_mode, wavelength_interval

ccd_gain = 3
ordinary_ray = 100
extra_ordinary_ray = 100
seeing = 1.5  # pix
_SPARC4_PLATE_SCALE = 0.35  # pix/arcsec
gaussian_std = seeing / _SPARC4_PLATE_SCALE
_SPARC4_POL_SEPARATION = 40
em_gain = ccd_operation_mode["em_gain"]
binn = ccd_operation_mode["binn"]
t_exp = ccd_operation_mode["t_exp"]
image_size = ccd_operation_mode["image_size"]
gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
shape = (image_size, image_size)
table = Table()
table["amplitude"] = [gaussian_amplitude]
table["x_mean"] = [image_size / 2]
table["y_mean"] = [image_size / 2]
table["x_stddev"] = [gaussian_std / binn]
table["y_stddev"] = [gaussian_std / binn]
table["theta"] = np.radians(np.array([0]))


@pytest.fixture
def chc1():
    return Concrete_Channel_1(
        wavelength_interval=wavelength_interval, ccd_operation_mode=ccd_operation_mode
    )


@pytest.fixture
def psf(chc1):
    return Point_Spread_Function(ccd_operation_mode, ccd_gain, seeing)


# ------------------------ Initialize the class --------------------------


def test_em_gain(psf):
    assert psf.em_gain == ccd_operation_mode["em_gain"]


def test_bin(psf):
    assert psf.binn == ccd_operation_mode["binn"]


def test_t_exp(psf):
    assert psf.t_exp == ccd_operation_mode["t_exp"]


def test_image_size(psf):
    assert psf.image_size == ccd_operation_mode["image_size"]


def test_ccd_gain(psf):
    assert psf.ccd_gain == ccd_gain


def test_make_noise_image(psf):
    new_noise_image = psf._make_noise_image(shape, table, gaussian_amplitude)
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) - gaussian_amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)

    assert np.allclose(
        noise_image, new_noise_image, atol=5 * np.sqrt(gaussian_amplitude)
    )


def test_calculate_star_psf_photometric(psf):
    star_coord = [image_size / 2, image_size / 2]
    table["amplitude"] = [gaussian_amplitude]
    star_image = make_gaussian_sources_image(shape, table)
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) - gaussian_amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)
    star_image += noise_image

    new_star_image = psf.create_star_psf(star_coord, ordinary_ray)

    assert np.allclose(star_image, new_star_image, atol=5 * np.sqrt(gaussian_amplitude))


def test_calculate_star_psf_photometric_polarimetric(psf):

    gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
    shape = (image_size, image_size)

    star_image = make_gaussian_sources_image(shape, table)
    star_coord = [
        image_size / 2 - _SPARC4_POL_SEPARATION,
        image_size / 2 - _SPARC4_POL_SEPARATION,
    ]
    gaussian_amplitude = extra_ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
    table["amplitude"] = [gaussian_amplitude]
    table["x_mean"] = [star_coord[0]]
    table["y_mean"] = [star_coord[1]]
    star_image += make_gaussian_sources_image(shape, table)
    psf_image = psf.create_star_psf(star_coord, ordinary_ray, extra_ordinary_ray)

    assert np.allclose(star_image, psf_image, atol=5 * np.sqrt(gaussian_amplitude))
