# -*- coding: utf-8 -*-
"""
Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import numpy as np
import pytest
from AIS.Atmosphere_Spectral_Response import Atmosphere_Spectral_Response

from .SPARC4_SR_curves import (
    specific_flux,
    wavelength_interval,
    wavelength_interval_len,
)


@pytest.fixture
def atm_sr():
    return Atmosphere_Spectral_Response()


# ---------------------------------------------------------------------------------------------------------


def test_read_spreadsheet(atm_sr):
    wavelength_interval = range(350, 1150, 50)
    transmitance = np.ones((1, wavelength_interval_len + 1))  # por hora vale
    atm_sr._read_spreadsheet()
    assert np.allclose(atm_sr.atm_wavelength_interval, wavelength_interval)
    assert np.allclose(atm_sr.transmitance, transmitance)


def test_calculate_spline(atm_sr):
    transmitance = np.ones((1, wavelength_interval_len))
    atm_sr._read_spreadsheet()
    new_transmitance = atm_sr._calculate_spline(wavelength_interval)
    assert np.allclose(new_transmitance, transmitance)


def test_apply_atmosphere_spectral_response(atm_sr):
    new_flux = atm_sr.apply_atmosphere_spectral_response(
        specific_flux, wavelength_interval
    )
    assert np.allclose(new_flux, specific_flux)
