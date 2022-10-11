# -*- coding: utf-8 -*-

"""
Flux Calculation class tests.

This script tests the operation of the Flux Calculation Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import numpy as np
import pytest
from AIS.Point_Spread_Function import Point_Spread_Function
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image, make_noise_image

from .AIS_spectral_response_curves import ccd_operation_mode


# _SPARC4_POL_SEPARATION = 40


@pytest.fixture
def psf():
    return Point_Spread_Function(ccd_operation_mode, channel=1)


# ------------------------ Initialize the class --------------------------


def test_ccd_operation_mode(psf):
    assert psf.ccd_operation_mode == ccd_operation_mode


def test_channel(psf):
    assert psf.channel == 1


def test_get_ccd_gain(psf):
    assert psf.ccd_gain == 3.37


# ------------------------------------------------------------------------------------------------------
seeing = 1
_SPARC4_PLATE_SCALE = 0.35
image_size = ccd_operation_mode['image_size']
half_img_size = image_size//2
shape = (image_size, image_size)
gaussian_std = seeing / _SPARC4_PLATE_SCALE
table = Table()
binn = ccd_operation_mode['binn']
x_coord = half_img_size
y_coord = half_img_size
table["x_mean"] = [x_coord]
table["y_mean"] = [y_coord]
table["x_stddev"] = [gaussian_std / binn]
table["y_stddev"] = [gaussian_std / binn]
table["theta"] = np.radians(np.array([0]))


def test_create_table(psf):
    psf._create_table(shape, seeing)
    assert psf.table == table


# ---------------------------------------------------------------------------------------------------------
ordinary_ray = 100
t_exp = ccd_operation_mode['t_exp']
em_gain = ccd_operation_mode['em_gain']
ccd_gain = 3.37
gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
table["amplitude"] = gaussian_amplitude


def test_make_noise_image(psf):
    image_noise = psf._make_noise_image(table, shape)
    new_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) -
        gaussian_amplitude
    )
    table["amplitude"] = [1]
    new_image *= make_gaussian_sources_image(shape, table)
    assert np.allclose(new_image, image_noise, atol=5 *
                       np.sqrt(gaussian_amplitude))


# ---------------------------------------------------------------------------------------------------------

star_coordinates = (half_img_size, half_img_size)
tmp_image_1 = make_gaussian_sources_image(shape, table)
amplitude = table["amplitude"]
table["amplitude"] = [1]
noise_image = (
    make_noise_image(shape, "poisson", amplitude) - amplitude
)
noise_image *= make_gaussian_sources_image(shape, table)
tmp_image_1 += noise_image


def test_create_image_ordinary_ray(psf):
    psf._create_table(star_coordinates, seeing)
    noise_image = psf._create_image_ordinary_ray(ordinary_ray)
    assert np.allclose(noise_image, tmp_image_1, atol=5 *
                       np.sqrt(gaussian_amplitude))

# ---------------------------------------------------------------------------------------------------------


_SPARC4_POL_SEPARATION = 40  # pix
table["x_mean"] -= _SPARC4_POL_SEPARATION
table["y_mean"] -= _SPARC4_POL_SEPARATION
table["amplitude"] = gaussian_amplitude
tmp_image_2 = make_gaussian_sources_image(shape, table)
table["amplitude"] = [1]
noise_image = (
    make_noise_image(shape, "poisson", gaussian_amplitude) - gaussian_amplitude
)
noise_image *= make_gaussian_sources_image(shape, table)
tmp_image_2 += noise_image
star_coordinates = (half_img_size, half_img_size)


def test_create_image_extra_ordinary_ray(psf):
    psf._create_table(star_coordinates, seeing)
    noise_image = psf._create_image_extra_ordinary_ray(ordinary_ray)
    assert np.allclose(noise_image, tmp_image_2, atol=5 *
                       np.sqrt(gaussian_amplitude))

# ---------------------------------------------------------------------------------------------------------


def test_creat_star_image(psf):
    star_image = psf.create_star_image(
        star_coordinates, ordinary_ray, ordinary_ray)
    assert np.allclose(star_image, tmp_image_2 + tmp_image_1, atol=5 *
                       np.sqrt(gaussian_amplitude))
