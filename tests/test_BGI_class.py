# -*- coding: utf-8 -*-
"""Tests of the Backgorund Image class.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""


import numpy as np
import pytest
from AIS.Background_Image import Background_Image

ccd_operation_mode = {
    "em_mode": 'Conv',
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 1024,
    'temp': -70
}


@pytest.fixture
def bgi_conv():
    return Background_Image(
        ccd_operation_mode, 1, bias_level=100
    )


@pytest.fixture
def bgi_em():
    dic = ccd_operation_mode.copy()
    dic["em_mode"] = 'EM'
    dic["em_gain"] = 2
    return Background_Image(dic, 1)


# ------------------------ Initialize the class --------------------------

def test_ccd_operation_mode(bgi_conv):
    assert bgi_conv.ccd_operation_mode == ccd_operation_mode


def test_bias_level(bgi_conv):
    assert bgi_conv.bias_level == 100


def test_channel(bgi_conv):
    assert bgi_conv.channel == 1


def test_noise_factor(bgi_conv):
    assert bgi_conv._NOISE_FACTOR == 1


def test_read_noise(bgi_conv):
    assert bgi_conv.read_noise == 6.67


def test_dark_noise(bgi_conv):
    assert bgi_conv.dark_noise == 5.8597559895090484e-05 * \
        ccd_operation_mode['t_exp']

# ----------------------------------------------------------------------------


def test_get_ccd_gain(bgi_conv):
    idx_tab = bgi_conv.get_ccd_gain()
    assert bgi_conv.ccd_gain == 3.37

# ----------------------- Calculate Background Image -------------------------


ccd_gain = 3.37
read_noise = 6.67
bias_level = 100


def test_create_bias_background(bgi_conv):
    bbg = bgi_conv.create_bias_background()
    bg_level = np.mean(bbg)
    new_noise = np.std(bbg)
    noise = read_noise / ccd_gain
    assert np.allclose(bg_level, bias_level, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)

# -------------------------------------------------------------------------------------------------------


dark_noise = 5.8597559895090484e-05 * ccd_operation_mode['t_exp']


def test_create_dark_background(bgi_conv):
    image = bgi_conv.create_dark_background()
    bg_level = np.mean(image)
    new_noise = np.std(image)
    noise = read_noise / ccd_gain
    assert np.allclose(bg_level, bias_level +
                       dark_noise / ccd_gain, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)

# -------------------------------------------------------------------------------------------------------


noise_factor = 1
_FLAT_LEVEL = 2**14
poisson_noise = _FLAT_LEVEL / ccd_gain
_PIXEL_SENSIBILITY = 0.03
em_gain = ccd_operation_mode['em_gain'] ** 2
binn = ccd_operation_mode['binn'] ** 2
noise = (
    np.sqrt(
        read_noise ** 2
        + (poisson_noise*(1+_PIXEL_SENSIBILITY))
        * noise_factor ** 2
        * em_gain ** 2
        * binn ** 2
    )
    / ccd_gain
)


def test_create_flat_image(bgi_conv):
    image = bgi_conv.create_flat_background()
    bg_level = np.mean(image)
    new_noise = np.std(image)
    assert np.allclose(bg_level, _FLAT_LEVEL, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)

# -------------------------------------------------------------------------------------------------------


sky_flux = 10
t_exp = ccd_operation_mode['t_exp']
background_level = bias_level + \
    (sky_flux * t_exp + dark_noise) * \
    em_gain * binn ** 2 / ccd_gain
noise = (
    np.sqrt(read_noise ** 2 + (sky_flux * t_exp + dark_noise) *
            noise_factor ** 2 * em_gain ** 2 * binn ** 2)
    / ccd_gain
)


def test_create_sky_background(bgi_conv):
    image = bgi_conv.create_sky_background(sky_flux)
    bg_level = np.mean(image)
    new_noise = np.std(image)
    assert np.allclose(bg_level, background_level, rtol=0.005)
    assert np.allclose(noise, new_noise, rtol=0.005)
