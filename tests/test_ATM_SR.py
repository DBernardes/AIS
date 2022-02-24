# -*- coding: utf-8 -*-
"""
Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import os

import numpy as np
import pandas as pd
import pytest
from AIS.Atmosphere_Spectral_Response import Atmosphere_Spectral_Response
from scipy.interpolate import splev, splrep

from .AIS_spectral_response_curves import (
    air_mass,
    atm_wavelength_interval,
    good_extinction_coef,
    photometric_extinction_coef,
    regular_extinction_coef,
    sky_condition,
    star_specific_photons_per_second,
    wavelength_interval,
)

star_specific_photons_per_second = star_specific_photons_per_second[0]


@pytest.fixture
def atm_sr():
    return Atmosphere_Spectral_Response(air_mass, sky_condition)


# ---------------------------------------------------------------------------------------------------------


def test_read_photometric_extinction_coef(atm_sr):
    atm_sr._read_spreadsheet()
    assert np.allclose(atm_sr.atm_wavelength_interval, atm_wavelength_interval)
    assert np.allclose(atm_sr.extinction_coefs, photometric_extinction_coef)


def test_read_regular_extinction_coef(atm_sr):
    atm_sr = Atmosphere_Spectral_Response(1, "regular")
    atm_sr._read_spreadsheet()
    assert np.allclose(atm_sr.atm_wavelength_interval, atm_wavelength_interval)
    assert np.allclose(atm_sr.extinction_coefs, regular_extinction_coef)


def test_read_good_extinction_coef(atm_sr):
    atm_sr = Atmosphere_Spectral_Response(1, "good")
    atm_sr._read_spreadsheet()
    assert np.allclose(atm_sr.atm_wavelength_interval, atm_wavelength_interval)
    assert np.allclose(atm_sr.extinction_coefs, good_extinction_coef)


def test_calculate_spline(atm_sr):
    transmitance = [10 ** (-0.4 * k * air_mass) for k in photometric_extinction_coef]
    spl = splrep(atm_wavelength_interval, transmitance)
    transmitance = splev(wavelength_interval, spl)
    atm_sr._read_spreadsheet()
    new_transmitance = atm_sr._calculate_atmosphere_transmitance(wavelength_interval)
    assert np.allclose(new_transmitance, transmitance)


def test_apply_atmosphere_spectral_response(atm_sr):
    atm_specific_photons_per_second = atm_sr.apply_atmosphere_spectral_response(
        star_specific_photons_per_second, wavelength_interval
    )

    transmitance = [10 ** (-0.4 * k * air_mass) for k in photometric_extinction_coef]
    spl = splrep(atm_wavelength_interval, transmitance)
    transmitance = splev(wavelength_interval, spl)
    star_specific_photons_per_second[0] = np.multiply(
        star_specific_photons_per_second[0], transmitance
    )

    assert np.allclose(
        atm_specific_photons_per_second, star_specific_photons_per_second
    )
