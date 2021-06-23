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


@pytest.fixture
def tel_sr():
    return Telescope_Spectral_Response()


l_init, l_final, l_step = 350, 1150, 50
wavelength_interv = range(l_init, l_final, l_step)
n = len(wavelength_interv)
specific_flux = np.ones((4, n))

ss = np.asarray(
    pd.read_excel(
        os.path.join("Telescope_Spectral_Response", "telescope_spectral_response.xlsx")
    )
)
tel_wavelength_interval = [float(value) for value in ss[1:, 0]]
reflectance = [float(value) for value in ss[1:, 1]]
tel_wavelength_interval = np.asarray(tel_wavelength_interval)
reflectance = np.asarray(reflectance)


# ---------------------------------------------------------------------------------------------------------------------


def test_read_spreadsheet(tel_sr):
    tel_sr._read_spreadsheet()
    assert np.allclose(tel_sr.tel_wavelength_interval, tel_wavelength_interval)
    assert np.allclose(tel_sr.reflectance, reflectance)


# def test_calculate_spline(tel_sr):
#     tel_sr._read_spreadsheet()
#     new_reflectance = tel_sr._calculate_spline(wavelength_interv)
#     spl = splrep(tel_wavelength_interval, reflectance)
#     reflectance = splev(wavelength_interval, spl)
#     assert np.allclose(new_reflectance, reflectance)


# def test_apply_telescope_spectral_response(tel_sr):
#     new_specific_flux = tel_sr.apply_telescope_spectral_response(
#         specific_flux, l_init, l_final, l_step
#     )
#     specific_flux = np.multiply(specific_flux, reflectance)
#     assert np.allclose(specific_flux, new_specific_flux)
