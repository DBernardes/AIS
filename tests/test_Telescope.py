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
from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength
from scipy.interpolate import splev, splrep
import pandas as pd

csv_file_name = 'telescope.csv'
BASE_PATH = os.path.join('AIS', 'Spectral_Response')
ss = pd.read_csv(os.path.join(BASE_PATH, csv_file_name))
wavelength = ss['Wavelength (nm)']
spectral_response = ss['Transmitance (%)']


@pytest.fixture
def telescope():
    return Telescope()


def test_csv_file_name(telescope):
    assert telescope.csv_file_name == csv_file_name
