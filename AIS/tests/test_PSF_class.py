# -*- coding: utf-8 -*-

"""
Flux Calculation class tests.

This script tests the operation of the Flux Calculation Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import numpy as np
import pytest
from AIS.Channel_Creator import Concrete_Channel_1
from AIS.Point_Spread_Function import Point_Spread_Function
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image

from .SPARC4_SR_curves import wavelength_interval

dic = {
    "em_gain": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 200,
}
ccd_gain = 3
ordinary_ray = 100
extra_ordinary_ray = 0
image_size = 200
star_coord = [100, 100]
seeing = 1.5  # pix
_SPARC4_PLATE_SCALE = 0.35  # pix/arcsec


@pytest.fixture
def chc1():
    return Concrete_Channel_1(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


@pytest.fixture
def psf(chc1):
    return Point_Spread_Function(chc1, dic, ccd_gain, seeing)


# ------------------------ Initialize the class --------------------------


def test_em_gain(psf):
    assert psf.em_gain == 1


def test_bin(psf):
    assert psf.binn == 1


def test_t_exp(psf):
    assert psf.t_exp == 1


def test_image_size(psf):
    assert psf.image_size == 200


def test_ccd_gain(psf):
    assert psf.ccd_gain == 3


# ----------------------- Calculate star psf -----------------------------


def test_calculate_star_psf_photometric(psf):
    em_gain = dic["em_gain"]
    binn = dic["binn"]
    t_exp = dic["t_exp"]
    gaussian_std = seeing * _SPARC4_PLATE_SCALE

    gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
    shape = (image_size, image_size)
    table = Table()
    table["amplitude"] = [gaussian_amplitude]
    table["x_mean"] = [star_coord[0]]
    table["y_mean"] = [star_coord[1]]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))

    star_image = make_gaussian_sources_image(shape, table)
    psf_image = psf.create_star_psf(star_coord, ordinary_ray)

    assert np.allclose(star_image, psf_image)


def test_calculate_star_psf_photometric_polarimetric(psf):
    em_gain = dic["em_gain"]
    binn = dic["binn"]
    t_exp = dic["t_exp"]
    gaussian_std = seeing * _SPARC4_PLATE_SCALE

    gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
    shape = (image_size, image_size)
    table = Table()
    table["amplitude"] = [gaussian_amplitude]
    table["x_mean"] = [star_coord[0]]
    table["y_mean"] = [star_coord[1]]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))
    star_image = make_gaussian_sources_image(shape, table)

    gaussian_amplitude = extra_ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
    table["amplitude"] = [gaussian_amplitude]
    table["x_mean"] = [star_coord[0]]
    table["y_mean"] = [star_coord[1]]
    star_image += make_gaussian_sources_image(shape, table)
    psf_image = psf.create_star_psf(star_coord, ordinary_ray, extra_ordinary_ray)

    assert np.allclose(star_image, psf_image)
