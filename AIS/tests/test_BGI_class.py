# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""


import numpy as np
import pytest
from AIS.Background_Image import Background_Image

ccd_operation_mode = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 1024,
}

em_gain = ccd_operation_mode["em_gain"]
binn = ccd_operation_mode["binn"]
t_exp = ccd_operation_mode["t_exp"]
preamp = ccd_operation_mode["preamp"]
hss = ccd_operation_mode["hss"]
image_size = ccd_operation_mode["image_size"]
sky_flux = 10
ccd_gain = 3
dark_current = 1e-5
read_noise = 6.67
bias_level = 500
dc = dark_current * t_exp
rn = read_noise
nf = 1
background_level = bias_level + (dc + sky_flux) * t_exp * em_gain * binn ** 2 / ccd_gain
noise = (
    np.sqrt(rn ** 2 + (sky_flux + dc) * t_exp * nf ** 2 * em_gain ** 2 * binn ** 2)
    / ccd_gain
)


@pytest.fixture
def bgi_conv():
    return Background_Image(
        ccd_operation_mode, ccd_gain, dark_current, read_noise, bias_level
    )


@pytest.fixture
def bgi_em():
    dic = ccd_operation_mode.copy()
    dic["em_mode"] = 1
    return Background_Image(dic, ccd_gain, dark_current, read_noise, bias_level)


# ------------------------ Initialize the class --------------------------


def test_em_gain(bgi_conv):
    assert bgi_conv.em_gain == em_gain


def test_preamp(bgi_conv):
    assert bgi_conv.preamp == preamp


def test_hss(bgi_conv):
    assert bgi_conv.hss == hss


def test_bin(bgi_conv):
    assert bgi_conv.binn == binn


def test_t_exp(bgi_conv):
    assert bgi_conv.t_exp == t_exp


def test_image_size(bgi_conv):
    assert bgi_conv.image_size == image_size


def test_ccd_gain(bgi_conv):
    assert bgi_conv.ccd_gain == ccd_gain


def test_bias_level(bgi_conv):
    assert bgi_conv.bias_level == bias_level


def test_noise_factor_1(bgi_conv):
    assert bgi_conv._NOISE_FACTOR == 1


def test_noise_factor_2(bgi_em):
    assert bgi_em._NOISE_FACTOR == 1.4


# ----------------------- Calculate Background Image -------------------------


def test_create_background_image(bgi_conv):
    image = bgi_conv.create_background_image(sky_flux)
    bg_level = np.mean(image)
    new_noise = np.std(image)
    assert np.allclose(bg_level, background_level, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)


def test_create_bias_image(bgi_conv):
    image = bgi_conv.create_bias_image()
    bg_level = np.mean(image)
    new_noise = np.std(image)
    noise = rn / ccd_gain
    assert np.allclose(bg_level, bias_level, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)


def test_create_dark_image(bgi_conv):
    image = bgi_conv.create_dark_image()
    bg_level = np.mean(image)
    new_noise = np.std(image)
    noise = rn / ccd_gain
    assert np.allclose(bg_level, bias_level + dc / ccd_gain, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)


def test_create_flat_image(bgi_conv):
    image = bgi_conv.create_flat_image()
    bg_level = np.mean(image)
    BIAS = 32000
    new_noise = np.std(image)
    poisson_noise = BIAS / ccd_gain
    pixel_sensibility_noise = BIAS / ccd_gain * 0.03  # 3% of pixel sensibility
    noise = (
        np.sqrt(
            rn ** 2
            + (poisson_noise + pixel_sensibility_noise)
            * nf ** 2
            * em_gain ** 2
            * binn ** 2
        )
        / ccd_gain
    )
    assert np.allclose(bg_level, BIAS, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)
