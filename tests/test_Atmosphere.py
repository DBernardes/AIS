# -*- coding: utf-8 -*-
"""
Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""


from AIS.Spectral_Response import Atmosphere
import pytest
import os
import numpy as np

from scipy.optimize import curve_fit
from scipy.interpolate import splev, splrep
import pandas as pd

obj_wavelength = np.linspace(400, 1100, 100)
csv_file_name = 'atmosphere.csv'
BASE_PATH = os.path.join('AIS', 'Spectral_Response')
spreadsheet_path = os.path.join(BASE_PATH, csv_file_name)
ss = pd.read_csv(spreadsheet_path)
wavelength = ss['Wavelength (nm)']
extinction_coef_photometric = ss['photometric']
extinction_coef_regular = ss['regular']
extinction_coef_good = ss['good']


def func(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d


@pytest.fixture
def atm():
    return Atmosphere()


def test_csv_file_name(atm):
    assert atm._CSV_FILE_NAME == csv_file_name


def test_base_path(atm):
    assert atm._BASE_PATH == BASE_PATH

# --------------------------------------------------------------------------------------------


def test_read_csv_file_photometric(atm):
    new_wavelength, new_transmitance = atm._read_csv_file(
        spreadsheet_path, 'photometric')
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, extinction_coef_photometric)


def test_read_csv_file_regular(atm):
    new_wavelength, new_transmitance = atm._read_csv_file(
        spreadsheet_path, 'regular')
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, extinction_coef_regular)


def test_read_csv_file_good(atm):
    new_wavelength, new_transmitance = atm._read_csv_file(
        spreadsheet_path, 'good')
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, extinction_coef_good)


air_mass = 1
spectral_response = 10 ** (-0.4 * extinction_coef_photometric * air_mass)
# spl = splrep(wavelength, spectral_response)
# new_spectral_response = splev(obj_wavelength, spl)
popt, _ = curve_fit(func, wavelength,
                    spectral_response)
new_spectral_response = func(obj_wavelength, *popt)


def test_interpolate_spectral_response(atm):
    class_spectral_response = atm._interpolate_spectral_response(
        wavelength, spectral_response, obj_wavelength)
    assert np.allclose(class_spectral_response, new_spectral_response)


sky_condition = 'photometric'


def test_get_spectral_response(atm):
    class_spectral_response = atm.get_spectral_response(
        obj_wavelength, air_mass, sky_condition)
    assert np.allclose(class_spectral_response, new_spectral_response)


esd = np.ones(len(obj_wavelength))
reduced_esd = np.multiply(new_spectral_response, esd)


def test_apply_spectral_response(atm):
    class_reduced_esd = atm.apply_spectral_response(
        esd, obj_wavelength, air_mass, sky_condition)
    assert np.allclose(class_reduced_esd, reduced_esd)
