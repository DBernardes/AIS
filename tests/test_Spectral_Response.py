# -*- coding: utf-8 -*-
"""
Tests of the Spectral Response Class

Oct 17th 2022

@author: denis
"""

from AIS.Spectral_Response import Spectral_Response
import pytest
import os
import numpy as np
from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength
from scipy.interpolate import splev, splrep
import pandas as pd

csv_file_name = 'csv_file.csv'
BASE_PATH = os.path.join('AIS', 'Spectral_Response')
ss = pd.read_csv(os.path.join(BASE_PATH, csv_file_name))
wavelength = ss['Wavelength (nm)']
spectral_response = ss['Transmitance (%)']/100


@pytest.fixture
def sr():
    return Spectral_Response()


def test_csv_file_name(sr):
    assert sr._CSV_FILE_NAME == csv_file_name


def test_base_path(sr):
    assert sr._BASE_PATH == BASE_PATH

# --------------------------------------------------------------------------------------------


def test_read_csv_file(sr):
    csv = os.path.join(BASE_PATH, csv_file_name)
    new_wavelength, new_transmitance = sr._read_csv_file(csv)
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, spectral_response)


spl = splrep(wavelength, spectral_response)
new_spectral_response = splev(obj_wavelength, spl)


def test_interpolate_spectral_response(sr):
    class_spectral_response = sr._interpolate_spectral_response(
        wavelength, spectral_response, obj_wavelength)
    assert np.allclose(class_spectral_response, new_spectral_response)


def test_get_spectral_response(sr):
    class_spectral_response = sr.get_spectral_response(obj_wavelength)
    assert np.allclose(class_spectral_response, new_spectral_response)


esd = np.ones(len(obj_wavelength))
reduced_esd = np.multiply(new_spectral_response, esd)


def test_apply_spectral_response(sr):
    class_reduced_esd = sr.apply_spectral_response(esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)
