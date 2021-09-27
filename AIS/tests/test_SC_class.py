# -*- coding: utf-8 -*-
"""Spectrum calculation tests.

This script tests the operation of the spectrum calculation class.
"""

import numpy as np
import pytest
from AIS.Spectrum_Calculation import Spectrum_Calculation

temperature = 5700
l_init = 400
l_final = 1150
l_step = 50


@pytest.fixture
def sc():
    return Spectrum_Calculation(temperature, l_init, l_final, l_step)


_H = 6.62607004e-34  # m2 kg / s
_C = 3e8  # m/s
_K = 1.38064852e-23  # m2 kg s-2 K-1
_telescope_effective_area = 0.804  # m2
_angular_aperture = 1

T = temperature
h = _H
c = _C
k = _K
specific_flux = []
for Lambda in range(l_init, l_final, l_step):
    Lambda *= 1e-9

    var1 = 2 * h * c ** 2 / Lambda ** 5
    var2 = np.e ** (h * c / (Lambda * k * T)) - 1
    black_body = var1 / var2
    photon_energy = h * c / Lambda
    photons_per_second = (
        black_body * _telescope_effective_area * _angular_aperture / photon_energy
    )

    specific_flux.append(photons_per_second)

specific_flux_length = len(specific_flux)
star_specific_flux = np.zeros((4, specific_flux_length))
star_specific_flux[0, :] = specific_flux


# ---------------------------------------Initialize the class ------------------------------------------------


def test_temperature(sc):
    assert sc.temperature == 5700


def test_l_init(sc):
    assert sc.l_init == 400


def test_l_final(sc):
    assert sc.l_final == 1150


def test_l_step(sc):
    assert sc.l_step == 50


# --------------------------------------------------------------------------------------------------------------


def test_calculate_sky_specific_flux(sc):
    sky_specific_flux = sc.calculate_sky_specific_flux()
    assert np.allclose(sky_specific_flux, star_specific_flux * 0.1)


def test_calculate_star_specific_flux(sc):
    star_specific_flux = sc.calculate_star_specific_flux()
    assert np.allclose(star_specific_flux, star_specific_flux)
