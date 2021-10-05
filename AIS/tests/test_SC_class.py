# -*- coding: utf-8 -*-
"""Spectrum calculation tests.

This script tests the operation of the spectrum calculation class.
"""

import numpy as np
import pytest
from AIS.Spectrum_Calculation import Spectrum_Calculation

from .SPARC4_SR_curves import wavelength_interval

temperature = 5700
l_init = 400
l_final = 1150
l_step = 50


@pytest.fixture
def sc():
    return Spectrum_Calculation(wavelength_interval, temperature)


_H = 6.62607004e-34  # m2 kg / s
_C = 3e8  # m/s
_K = 1.38064852e-23  # m2 kg s-2 K-1
B = 0.2e-6  # m
S_0 = 4e-2  # W/m2/m
tel_area = 0.804  # m2
magnitude = 22

T = temperature
h = _H
c = _C
k = _K
temp = []
num = int((l_final - l_init) / l_step)
for Lambda in np.linspace(l_init, l_final, num):
    Lambda *= 1e-9
    photons_number = S_0 * 10 ** (-magnitude / 2.5) * Lambda * B * tel_area / (h * c)
    temp.append(photons_number)

specific_flux = np.zeros((4, num))
specific_flux[0, :] = temp

# ---------------------------------------Initialize the class ------------------------------------------------


def test_temperature(sc):
    assert sc.star_temperature == 5700


def test_l_init(sc):
    assert sc.wavelength_interval == wavelength_interval


# --------------------------------------------------------------------------------------------------------------
def test_calculate_specific_flux(sc):
    star_specific_flux = sc.calculate_specific_flux(magnitude)
    assert np.allclose(star_specific_flux, star_specific_flux)
