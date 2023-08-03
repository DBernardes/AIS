# -*- coding: utf-8 -*-
"""
Tests of the Spectral Response Class

Oct 19th 2022

@author: denis
"""

from AIS.Spectral_Response import Telescope
import pytest
import os
import numpy as np
from scipy.interpolate import splev, splrep, PchipInterpolator
import pandas as pd
from copy import copy

obj_wavelength = np.linspace(400, 1100, 100)
csv_file_name = "telescope.csv"
BASE_PATH = os.path.join("AIS", "Spectral_Response")
spreadsheet_path = os.path.join(BASE_PATH, csv_file_name)
ss = pd.read_csv(spreadsheet_path)
wavelength = ss["Wavelength (nm)"]
spectral_response = ss["Transmitance (%)"] / 100
_PRIMARY_MIRROR_ADJUSTMENT = 0.933
_SECONDARY_MIRROR_ADJUSTMENT = 0.996
spl = PchipInterpolator(wavelength, spectral_response)
new_spectral_response = (
    spl(obj_wavelength) ** 2 * _PRIMARY_MIRROR_ADJUSTMENT * _SECONDARY_MIRROR_ADJUSTMENT
)


@pytest.fixture
def telescope():
    return Telescope()


def test_csv_file_name(telescope):
    assert telescope._CSV_FILE_NAME == csv_file_name


def test_base_path(telescope):
    assert telescope._BASE_PATH == BASE_PATH


# --------------------------------------------------------------------------------------------


def test_read_csv_file(telescope):
    new_wavelength, new_transmitance = telescope._read_csv_file(spreadsheet_path)
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, spectral_response)


def test_interpolate_spectral_response(telescope):
    telescope.obj_wavelength = obj_wavelength
    class_spectral_response = telescope._interpolate_spectral_response(
        wavelength, spectral_response
    )
    assert np.allclose(class_spectral_response, spl(obj_wavelength), atol=1e-3)


def test_get_spectral_response(telescope):
    class_spectral_response = telescope.get_spectral_response(obj_wavelength)
    assert np.allclose(class_spectral_response, new_spectral_response)


sed = np.ones((4, len(obj_wavelength)))


def test_apply_spectral_response(telescope):
    class_reduced_esd = telescope.apply_spectral_response(obj_wavelength, copy(sed))
    new_sed = np.multiply(new_spectral_response, sed)
    assert np.allclose(class_reduced_esd, new_sed)
