# -*- coding: utf-8 -*-
"""
Tests of the Spectral Response Class

Oct 17th 2022

@author: denis
"""

from AIS.Spectral_Response import Spectral_Response
import pytest, unittest
import os
import numpy as np
from scipy.interpolate import splev, splrep, interp1d, PchipInterpolator
import pandas as pd
from copy import copy


class Test_Spectral_Response(unittest.TestCase):

    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.ones(100)
    BASE_PATH = os.path.join("AIS", "Spectral_Response")
    CSV_FILE_NAME = "csv_file.csv"

    @classmethod
    def setUpClass(cls):
        cls.sr = Spectral_Response()
        spreadsheet_path = os.path.join(cls.BASE_PATH, cls.CSV_FILE_NAME)
        ss = pd.read_csv(spreadsheet_path)
        wavelength = ss["Wavelength (nm)"]
        spectral_response = ss["Transmitance (%)"] / 100
        cls.wavelength, cls.spectral_response = wavelength, spectral_response
    
        spl = PchipInterpolator(wavelength, spectral_response)
        cls.interpolated_spectral_response = spl(cls.OBJ_WAVELENGTH) 

    def test_csv_file_name(self):
        assert self.sr.CSV_FILE_NAME == self.CSV_FILE_NAME

    def test_base_path(self):
        assert self.sr.BASE_PATH == self.BASE_PATH

    def test_read_csv_file(self):
        csv = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)
        new_wavelength, new_transmitance = self.sr._read_csv_file(csv)
        assert np.allclose(new_wavelength, self.wavelength)
        assert np.allclose(new_transmitance, self.spectral_response)

    def test_interpolate_spectral_response(self):
        sr = self.sr
        sr.obj_wavelength = self.OBJ_WAVELENGTH
        class_spectral_response = sr._interpolate_spectral_response(
            self.wavelength, self.spectral_response
        )
        assert np.allclose(class_spectral_response, self.interpolated_spectral_response)

    def test_get_spectral_response(self):
        new_spectral_response = self.sr.get_spectral_response(self.OBJ_WAVELENGTH)
        assert np.allclose(new_spectral_response, self.interpolated_spectral_response)

    def test_check_var_in_a_list(self):
        self.sr._check_var_in_a_list('dummy', 'dummy', ['dummy'])

    def test_check_var_in_a_list_err(self):
        with pytest.raises(ValueError):
            self.sr._check_var_in_a_list('dummy', 'dummy', [])

    def test_verify_var_interv(self):
        self.sr._verify_var_in_interval(10, 'dummy')
    
    def test_verify_var_interv_err(self):
        with pytest.raises(ValueError):
            self.sr._verify_var_in_interval(-1, 'dummy')

    def test_apply_spectral_response(self):
        spectral_response = self.sr.get_spectral_response(self.OBJ_WAVELENGTH)
        sed = np.multiply(spectral_response, self.SED)
        new_sed = self.sr.apply_spectral_response(self.OBJ_WAVELENGTH, self.SED)
        
        assert np.allclose(sed, new_sed)
