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

csv_file_name = 'csv_file.csv'
BASE_PATH = os.path.join('AIS', 'Spectral_Response')


@pytest.fixture
def sr():
    return Spectral_Response(csv_file_name)


def test_csv_file_name(sr):
    assert sr.csv_file_name == csv_file_name


def test_base_path(sr):
    assert sr._BASE_PATH == BASE_PATH

# --------------------------------------------------------------------------------------------


def test_read_csv_file(sr):
    csv = os.path.join(BASE_PATH, csv_file_name)
    wavelength, transmitance = sr._read_csv_file(csv)
    assert np.allclose(wavelength, [0])
    assert np.allclose(transmitance, [0])


def test_get_spectral_response(sr):
    assert np.allclose(sr.get_spectral_response(), [0])
