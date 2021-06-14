# -*- coding: utf-8 -*-
"""Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import numpy as np
import pytest
from AIS.Atmosphere_Spectral_Response import Atmosphere_Spectral_Response

l_init, l_final, l_step = 350, 1100, 50
wavelength_interv = np.asarray(range(l_init, l_final, l_step))
n = len(wavelength_interv)
specific_flux = np.ones((4, n))


@pytest.fixture
def atm_sr():
    return Atmosphere_Spectral_Response()


def test_read_spreadsheet(atm_sr):
    atm_sr._read_spreadsheet()
    assert atm_sr.atm_wavelength_interval.all() == wavelength_interv.all()
    assert atm_sr.transmitance.all() == specific_flux[0, :].all()


def test_calculate_spline(atm_sr):
    transmitance = np.ones((1, n))[0]
    atm_sr._read_spreadsheet()
    atm_sr._calculate_spline(wavelength_interv)
    assert atm_sr.transmitance.all() == transmitance.all()


def test_apply_atmosphere_spectral_response(atm_sr):
    new_flux = atm_sr.apply_atmosphere_spectral_response(
        specific_flux, l_init, l_final, l_step
    )
    assert new_flux.all() == specific_flux.all()
