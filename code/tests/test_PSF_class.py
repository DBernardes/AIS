# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Flux Calculation Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


from code.CHC.CHC import Concrete_Channel_1
from code.PSF.PSF import Point_Spread_Function

import numpy as np
import pytest
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 200,
}


@pytest.fixture
def chc1():
    return Concrete_Channel_1(ccd_temp=-70, sparc4_operation_mode="phot")


@pytest.fixture
def psf(chc1):
    return Point_Spread_Function(chc1, dic, 3)


# ------------------------ Initialize the class --------------------------


def test_em_gain(psf):
    assert psf.em_gain == 1


def test_bin(psf):
    assert psf.binn == 1


def test_t_exp(psf):
    assert psf.t_exp == 1


def test_ccd_gain(psf):
    assert psf.ccd_gain == 3


# ----------------------- Calculate star PSF -----------------------------


def test_calculate_star_PSF(psf):
    em_gain = dic["em_gain"]
    binn = dic["binn"]
    t_exp = dic["t_exp"]
    ccd_gain = 3
    gaussian_std = 3
    star_flux = 100
    image_size = 200
    star_coord = [100, 100]

    gaussian_amplitude = star_flux * t_exp * em_gain * binn ** 2 / ccd_gain
    shape = (image_size, image_size)
    table = Table()
    table["amplitude"] = [gaussian_amplitude]
    table["x_mean"] = [star_coord[0]]
    table["y_mean"] = [star_coord[1]]
    table["x_stddev"] = [gaussian_std / binn]
    table["y_stddev"] = [gaussian_std / binn]
    table["theta"] = np.radians(np.array([0]))

    star_image = make_gaussian_sources_image(shape, table)

    assert np.sum(psf.create_star_PSF(star_flux, star_coord, gaussian_std)) == np.sum(
        star_image
    )
