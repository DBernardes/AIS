# -*- coding: utf-8 -*-

"""Test of the Telescope Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import os

import numpy as np
import pandas as pd
import pytest
from AIS.Telescope_Spectral_Response import Telescope_Spectral_Response
from scipy.interpolate import splev, splrep

from .AIS_spectral_response_curves import (
    photometric_extinction_coef,
    star_specific_photons_per_second,
    tel_reflectance,
    tel_reflectance_spreadsheet,
    tel_wavelength_interval,
    wavelength_interval,
    wavelength_interval_len,
)


@pytest.fixture
def tel_sr():
    return Telescope_Spectral_Response()


# ---------------------------------------------------------------------------------------------------------------------


def test_read_spreadsheet(tel_sr):
    tel_sr._read_spreadsheet()
    assert np.allclose(tel_sr.tel_wavelength_interval, tel_wavelength_interval)
    assert np.allclose(tel_sr.reflectance, tel_reflectance_spreadsheet)


def test_calculate_spline(tel_sr):
    ss = pd.read_csv(
        os.path.join(
            "AIS", "Telescope_Spectral_Response", "telescope_spectral_response.csv"
        ),
        dtype=np.float64,
        skiprows=1,
        decimal=".",
    )
    tel_wavelength_interval = ss["(nm)"]
    reflectance = ss["(%)"] / 100

    tel_sr._read_spreadsheet()
    new_reflectance = tel_sr._calculate_spline(wavelength_interval)
    spl = splrep(tel_wavelength_interval, reflectance)
    reflectance = splev(wavelength_interval, spl)
    assert np.allclose(new_reflectance, reflectance)


def test_apply_telescope_spectral_response(tel_sr):
    ss = pd.read_csv(
        os.path.join(
            "AIS", "Telescope_Spectral_Response", "telescope_spectral_response.csv"
        ),
        dtype=np.float64,
        skiprows=1,
        decimal=".",
    )
    tel_wavelength_interval = ss["(nm)"]
    reflectance = ss["(%)"] / 100

    spl = splrep(tel_wavelength_interval, reflectance)
    reflectance = splev(wavelength_interval, spl)
    new_star_specific_photons_per_second = np.multiply(
        star_specific_photons_per_second, reflectance
    )

    new_star_specific_photons_per_second = tel_sr.apply_telescope_spectral_response(
        star_specific_photons_per_second, wavelength_interval
    )
    assert np.allclose(
        star_specific_photons_per_second, new_star_specific_photons_per_second
    )