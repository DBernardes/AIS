# -*- coding: utf-8 -*-
"""
Tests of the Spectral Response Class

Oct 19th 2022

@author: denis
"""

import os
import unittest
from copy import copy

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import PchipInterpolator, splev, splrep

from AIS.Spectral_Response import Telescope


class Test_Telescope(unittest.TestCase):
    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.ones(100)
    ADJUSMENT_COEFFS = {"20230328": (0.86, 0.95), "20230329": (0.998, 0.995)}
    BASE_PATH = os.path.join("AIS", "Spectral_Response", "telescope")

    @classmethod
    def setUpClass(cls):
        cls.telescope = Telescope()
        spreadsheet_path = os.path.join(cls.BASE_PATH, "aluminium_curve.csv")
        ss = pd.read_csv(spreadsheet_path)
        wavelength = ss["Wavelength (nm)"]
        spectral_response = ss["Transmitance (%)"] / 100
        cls.wavelength, cls.spectral_response = wavelength, spectral_response

        spl = PchipInterpolator(wavelength, spectral_response)
        cls.interpolated_spectral_response = spl(cls.OBJ_WAVELENGTH)

    def test_csv_file_name(self):
        assert self.telescope.CSV_FILE_NAME == "aluminium_curve.csv"

    def test_base_path(self):
        assert self.telescope.BASE_PATH == self.BASE_PATH

    # --------------------------------------------------------------------------------------------

    def test_read_csv_file(self):
        spreadsheet_path = os.path.join(self.BASE_PATH, "aluminium_curve.csv")
        new_wavelength, new_transmitance = self.telescope._read_csv_file(
            spreadsheet_path
        )
        assert np.allclose(new_wavelength, self.wavelength)
        assert np.allclose(new_transmitance, self.spectral_response)

    def test_interpolate_spectral_response(self):
        telescope = self.telescope
        telescope.obj_wavelength = self.OBJ_WAVELENGTH
        obj_spectral_response = telescope._interpolate_spectral_response(
            self.wavelength, self.spectral_response
        )
        assert np.allclose(
            obj_spectral_response, self.interpolated_spectral_response, atol=1e-3
        )

    def test_get_spectral_response(self):
        primary_coeff, secondary_coeff = self.ADJUSMENT_COEFFS["20230329"]
        spectral_response = (
            self.interpolated_spectral_response**2 * primary_coeff * secondary_coeff
        )
        new_spectral_response = self.telescope.get_spectral_response(
            self.OBJ_WAVELENGTH
        )
        assert np.allclose(new_spectral_response, spectral_response)

    def test_apply_spectral_response(self):
        spectral_response = self.telescope.get_spectral_response(self.OBJ_WAVELENGTH)
        sed = np.multiply(spectral_response, self.SED)
        new_sed = self.telescope.apply_spectral_response(self.OBJ_WAVELENGTH, self.SED)
        assert np.allclose(new_sed, sed)
