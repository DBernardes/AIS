# -*- coding: utf-8 -*-
"""Spectrum calculation tests.

This script tests the operation of the spectrum calculation class.
"""

import numpy as np
import pytest
from AIS.Spectrum_Calculation import Spectrum_Calculation

from .AIS_spectral_response_curves import wavelength_interval, wavelength_interval_len

temperature = 5700


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
for Lambda in wavelength_interval:
    Lambda *= 1e-9
    photons_number = S_0 * 10 ** (-magnitude / 2.5) * Lambda * B * tel_area / (h * c)
    temp.append(photons_number)

specific_flux = np.zeros((4, wavelength_interval_len))
specific_flux[0, :] = temp

# ---------------------------------------Initialize the class ------------------------------------------------


def test_temperature(sc):
    assert sc.star_temperature == 5700


def test_l_init(sc):
    assert sc.wavelength_interval == wavelength_interval


# --------------------------------------------------------------------------------------------------------------
def test_calculate_specific_photons_per_second(sc):
    star_specific_photons_per_second = sc.calculate_specific_photons_per_second(
        magnitude
    )
    assert np.allclose(
        star_specific_photons_per_second, star_specific_photons_per_second
    )