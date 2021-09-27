# -*- coding: utf-8 -*-
"""Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import numpy as np
import pytest
from AIS.Atmosphere_Spectral_Response import Atmosphere_Spectral_Response

l_init, l_final, l_step = 400, 1150, 50
wavelength_interv = np.asarray(range(l_init, l_final, l_step))
n = len(wavelength_interv)
specific_flux = np.ones((4, n))


@pytest.fixture
def atm_sr():
    return Atmosphere_Spectral_Response()


# ---------------------------------------------------------------------------------------------------------


def test_read_spreadsheet(atm_sr):
    wavelength_interv = range(350, 1150, 50)  # por hora vale
    n = len(wavelength_interv)
    transmitance = np.ones((1, n))  # por hora vale
    atm_sr._read_spreadsheet()
    assert np.allclose(atm_sr.atm_wavelength_interval, wavelength_interv)
    assert np.allclose(atm_sr.transmitance, transmitance)


def test_calculate_spline(atm_sr):
    transmitance = np.ones((1, n))
    atm_sr._read_spreadsheet()
    new_transmitance = atm_sr._calculate_spline(wavelength_interv)
    assert np.allclose(new_transmitance, transmitance)


def test_apply_atmosphere_spectral_response(atm_sr):
    new_flux = atm_sr.apply_atmosphere_spectral_response(
        specific_flux, l_init, l_final, l_step
    )
    assert np.allclose(new_flux, specific_flux)
