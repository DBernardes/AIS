# -*- coding: utf-8 -*-
"""Spectrum calculation tests.

This script tests the operation of the spectrum calculation class.
"""

import numpy as np
import pytest
from AIS.Spectrum_Calculation import Spectrum_Calculation

temperature = 5700
l_init = 350
l_final = 1150
l_step = 50


@pytest.fixture
def sc():
    return Spectrum_Calculation(temperature, l_init, l_final, l_step)


def calculate_specific_flux(temperature, l_init, l_final, l_step):
    h = 6.62607004e-34  # m2 kg / s
    c = 3e8  # m/s
    k = 1.38064852e-23  # m2 kg s-2 K-1
    T = temperature
    specific_flux = []
    for Lambda in range(l_init, l_final, l_step):
        Lambda *= 1e-9
        B = 2 * h * c ** 2 / Lambda ** 5 * 1 / (np.e ** (h * c / (Lambda * k * T)) - 1)
        specific_flux.append(B)
    return np.asarray(specific_flux)


temp = calculate_specific_flux(temperature, l_init, l_final, l_step)
n = len(temp)
specific_flux = np.zeros((4, n))
specific_flux[0, :] = temp  # usar este

# ---------------------------------------Initialize the class ------------------------------------------------


def test_temperature(sc):
    assert sc.temperature == 5700


def test_l_init(sc):
    assert sc.l_init == 350


def test_l_final(sc):
    assert sc.l_final == 1150


def test_l_step(sc):
    assert sc.l_step == 50


# --------------------------------------------------------------------------------------------------------------


def test_calculate_sky_specific_flux(sc):
    sky_specific_flux = sc.calculate_sky_specific_flux()
    assert np.allclose(sky_specific_flux, specific_flux * 0.1)


def test_calculate_star_specific_flux(sc):
    star_specific_flux = sc.calculate_star_specific_flux()
    assert np.allclose(star_specific_flux, specific_flux)
