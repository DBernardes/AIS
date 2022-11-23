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
from math import pi

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
star_coordinates = (half_img_size, half_img_size)
gaussian_std = seeing / _SPARC4_PLATE_SCALE
binn = ccd_operation_mode['binn']
x_coord = half_img_size
y_coord = half_img_size


def test_create_table(psf):
    table = Table()
    table["x_mean"] = [x_coord]
    table["y_mean"] = [y_coord]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    psf._create_table(star_coordinates, seeing)
    assert psf.table == table


# ---------------------------------------------------------------------------------------------------------

shape = (image_size, image_size)
ordinary_ray = 100
t_exp = ccd_operation_mode['t_exp']
em_gain = ccd_operation_mode['em_gain']
ccd_gain = 3.37
gaussian_amplitude = ordinary_ray * t_exp * \
    em_gain * binn ** 2 / (ccd_gain * 2 * pi)


def test_make_noise_image(psf):
    table = Table()
    table["x_mean"] = [x_coord]
    table["y_mean"] = [y_coord]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    table["amplitude"] = gaussian_amplitude
    psf._create_table(star_coordinates, seeing)
    psf.table["amplitude"] = gaussian_amplitude
    image_noise = psf._make_noise_image()
    new_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) -
        gaussian_amplitude
    )
    table["amplitude"] = [1]
    new_image *= make_gaussian_sources_image(shape, table)
    assert np.allclose(new_image, image_noise, atol=5 *
                       np.sqrt(gaussian_amplitude))


# ---------------------------------------------------------------------------------------------------------


def test_create_image_ordinary_ray(psf):
    table = Table()
    table["x_mean"] = [x_coord]
    table["y_mean"] = [y_coord]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    table["amplitude"] = gaussian_amplitude
    tmp_image_1 = make_gaussian_sources_image(shape, table)
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) -
        gaussian_amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)
    psf._create_table(star_coordinates, seeing)
    star_image = psf._create_image_ordinary_ray(gaussian_amplitude)

    assert np.allclose(star_image, tmp_image_1 + noise_image, atol=5 *
                       np.sqrt(gaussian_amplitude))


# # ---------------------------------------------------------------------------------------------------------
_SPARC4_POL_SEPARATION = 20  # pix


def test_create_image_extra_ordinary_ray(psf):
    table = Table()
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    table["x_mean"] = [x_coord - _SPARC4_POL_SEPARATION]
    table["y_mean"] = [y_coord - _SPARC4_POL_SEPARATION]
    table["amplitude"] = gaussian_amplitude
    tmp_image_2 = make_gaussian_sources_image(shape, table)
    amplitude = table["amplitude"]
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", amplitude) -
        amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)
    tmp_image_2 += noise_image
    psf._create_table(star_coordinates, seeing)
    noise_image = psf._create_image_extra_ordinary_ray(ordinary_ray)
    assert np.allclose(noise_image, tmp_image_2, atol=5 *
                       np.sqrt(gaussian_amplitude))

# ---------------------------------------------------------------------------------------------------------


def test_creat_star_image(psf):
    table = Table()
    table["x_mean"] = [x_coord]
    table["y_mean"] = [y_coord]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    table["amplitude"] = gaussian_amplitude
    tmp_image_1 = make_gaussian_sources_image(shape, table)
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) -
        gaussian_amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)
    tmp_image_1 += noise_image

    table["x_mean"] = [x_coord - _SPARC4_POL_SEPARATION]
    table["y_mean"] = [y_coord - _SPARC4_POL_SEPARATION]
    table["amplitude"] = gaussian_amplitude
    tmp_image_2 = make_gaussian_sources_image(shape, table)
    table["amplitude"] = [1]
    noise_image = (
        make_noise_image(shape, "poisson", gaussian_amplitude) -
        gaussian_amplitude
    )
    noise_image *= make_gaussian_sources_image(shape, table)
    tmp_image_2 += noise_image

    star_image = psf.create_star_image(
        star_coordinates, ordinary_ray, ordinary_ray)
    assert np.allclose(star_image, tmp_image_1 + tmp_image_2, atol=5 *
                       np.sqrt(gaussian_amplitude))
