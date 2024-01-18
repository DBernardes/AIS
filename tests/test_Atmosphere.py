# -*- coding: utf-8 -*-
"""
Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""


from AIS.Spectral_Response import Atmosphere
import pytest, unittest
import os
import numpy as np
from copy import copy
from scipy.optimize import curve_fit
from scipy.interpolate import splev, splrep, interp1d, PchipInterpolator
import pandas as pd

class Test_Atmosphere(unittest.TestCase):

    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.ones(100)
    BASE_PATH = os.path.join("AIS", "Spectral_Response")
    CSV_FILE_NAME = "atmosphere_profile.csv"
    ATM_EXTINCTION_FILE = "willton.csv"
    AIR_MASS = 1
    SKY_CONDITION = "photometric"

    @classmethod
    def setUpClass(cls):
        cls.atm = Atmosphere()
        spreadsheet_path = os.path.join(cls.BASE_PATH, cls.CSV_FILE_NAME)
        ss = pd.read_csv(spreadsheet_path)
        wavelength = ss["Wavelength (nm)"]
        spectral_response = ss["Transmitance (%)"] / 100
        cls.wavelength, cls.spectral_response = wavelength, spectral_response
    
        spl = PchipInterpolator(wavelength, spectral_response)
        cls.interpolated_spectral_response = spl(cls.OBJ_WAVELENGTH) 



    def test_csv_file_name(self):
        assert self.atm.CSV_FILE_NAME == self.CSV_FILE_NAME


    def test_base_path(self):
        assert self.atm.BASE_PATH == self.BASE_PATH


    def test_read_csv_file(self):
        spreadsheet_path = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)
        new_wavelength, new_transmitance = self.atm._read_csv_file(spreadsheet_path)
        assert np.allclose(new_wavelength, self.wavelength)
        assert np.allclose(new_transmitance, self.spectral_response)


    def test_interpolate_spectral_response(self):
        atm = self.atm
        atm.obj_wavelength = self.OBJ_WAVELENGTH
        class_spectral_response = atm._interpolate_spectral_response(
            self.wavelength, self.spectral_response
        )
        assert np.allclose(class_spectral_response, self.interpolated_spectral_response)


    def test_atmosphere_profile_func(self):
        C = 2
        wavelength_interv = np.linspace(400, 1000, 50)
        csv_file = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)

        sys_wavelength, spectral_response = self.atm._read_csv_file(csv_file)
        spl = interp1d(sys_wavelength, spectral_response, bounds_error=False, kind="cubic")
        spectral_response = spl(wavelength_interv) * C

        new_spectral_response = self.atm._func(wavelength_interv, C)

        assert np.allclose(spectral_response, new_spectral_response)

    def test_get_spectral_response(self):
        extinction_ss = pd.read_csv(os.path.join(self.BASE_PATH, self.ATM_EXTINCTION_FILE))
        extinction_coef = 10 ** (-0.4 * self.AIR_MASS * extinction_ss[self.SKY_CONDITION])
        popt_M1, _ = curve_fit(self.atm._func, extinction_ss["Wavelength (nm)"], extinction_coef)
        new_spectral_response = self.interpolated_spectral_response * popt_M1[0]

        spectral_response = self.atm.get_spectral_response(
            self.OBJ_WAVELENGTH, self.AIR_MASS, self.SKY_CONDITION
        )
        
        

        assert np.allclose(new_spectral_response, spectral_response)


    def test_apply_spectral_response(self):
        spectral_response = self.atm.get_spectral_response(
            self.OBJ_WAVELENGTH, self.AIR_MASS, self.SKY_CONDITION
        )
        sed = np.multiply(spectral_response, self.SED)
        new_sed = self.atm.apply_spectral_response(
            self.OBJ_WAVELENGTH, self.SED, self.AIR_MASS, self.SKY_CONDITION
        )

        assert np.allclose(sed, new_sed)
