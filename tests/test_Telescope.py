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

class Test_Telescope:

    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    PRIMARY_MIRROR_ADJUSTMENT = 0.933
    SECONDARY_MIRROR_ADJUSTMENT = 0.996
    BASE_PATH = os.path.join("AIS", "Spectral_Response")
    telescope = Telescope()


    def __init__(self):
        spreadsheet_path = os.path.join(self.BASE_PATH, 'telescope.csv')
        ss = pd.read_csv(spreadsheet_path)
        wavelength = ss["Wavelength (nm)"]
        spectral_response = ss["Transmitance (%)"] / 100
    
        spl = PchipInterpolator(wavelength, spectral_response)
        self.new_spectral_response = (
            spl(self.OBJ_WAVELENGTH) ** 2 * self.PRIMARY_MIRROR_ADJUSTMENT * self.SECONDARY_MIRROR_ADJUSTMENT
        )


    def test_csv_file_name(self):
        assert self.telescope._CSV_FILE_NAME == "telescope.csv"


    def test_base_path(self):
        assert telescope._BASE_PATH == self.BASE_PATH


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
