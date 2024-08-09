# # -*- coding: utf-8 -*-
# Tests of the Channel Class
# Oct 17th 2022
# @author: denis
#

import os
import unittest
from copy import copy
from math import atan, sqrt

import numpy as np
import pandas as pd
import pytest
from numpy import cos, pi, sin
from scipy.interpolate import PchipInterpolator, interp1d, splev, splrep

from AIS.Spectral_Response import Channel
from AIS.Spectral_Response._utils import (
    apply_matrix,
    calculate_depolarizer_matrix,
    calculate_polarizer_matrix,
    calculate_retarder_matrix,
)


class Test_Channel(unittest.TestCase):
    CHANNEL = 1
    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.zeros((4, 100))
    SED[0, :] = 1
    BASE_PATH = os.path.join("AIS", "Spectral_Response")
    CSV_FILE_NAME = "polarizer.csv"
    BASE_PATH = os.path.join("AIS", "Spectral_Response", "channel")
    POL_OPTICAL_COMPONENTS = {
        "retarder": "retarder.csv",
        "analyzer": "analyzer.csv",
    }
    PHOT_OPTICAL_COMPONENTS = {
        "collimator": "collimator.csv",
        "dichroic": "dichroic.csv",
        "camera": "camera.csv",
        "ccd": "ccd.csv",
    }
    POLARIZER_ANGLE = 0

    @classmethod
    def setUpClass(cls):
        cls.channel = Channel(cls.CHANNEL)
        spreadsheet_path = os.path.join(cls.BASE_PATH, cls.CSV_FILE_NAME)
        ss = pd.read_csv(spreadsheet_path)
        wavelength = ss["Wavelength (nm)"]
        spectral_response = ss["Transmitance (%)"] / 100
        cls.wavelength, cls.spectral_response = wavelength, spectral_response

        spl = splrep(
            wavelength,
            spectral_response,
        )
        interpolated_spectral_response = splev(cls.OBJ_WAVELENGTH, spl)
        interpolated_spectral_response[interpolated_spectral_response < 0] = 0
        cls.interpolated_spectral_response = interpolated_spectral_response

    # def _apply_matrix(matrix, sed):
    #     for idx, _ in enumerate(sed[0]):
    #         sed[:, idx] = np.transpose(matrix.dot(sed[:, idx]))
    #     return sed

    def test_csv_file_name_phot(self):
        assert self.channel.PHOT_OPTICAL_COMPONENTS == self.PHOT_OPTICAL_COMPONENTS

    def test_csv_file_name_pol(self):
        assert self.channel.POL_OPTICAL_COMPONENTS == self.POL_OPTICAL_COMPONENTS

    # def test_base_path(self):
    #     assert self.channel.BASE_PATH == self.BASE_PATH

    def test_read_csv_file_polarizer(self):
        spreadsheet_path = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)

        new_wavelength, new_spectral_response = self.channel._read_csv_file(
            spreadsheet_path
        )
        assert np.allclose(new_wavelength, self.wavelength)
        assert np.allclose(new_spectral_response, self.spectral_response)

    def test_read_csv_file(self):
        components_files = [
            "collimator.csv",
            "analyzer.csv",
            "collimator.csv",
            "depolarizer.csv",
            "polarizer.csv",
            "polarizer_contrast_ratio.csv",
            "retarder.csv",
            "retarder_phase_diff_half.csv",
            "retarder_phase_diff_quarter.csv",
        ]
        for csv_file in components_files:
            spreadsheet_path = os.path.join(self.BASE_PATH, csv_file)
            pd.read_csv(spreadsheet_path)

    def test_read_csv_file_2(self):
        components_files = ["dichroic", "camera", "ccd"]
        for channel_id in [1, 2, 3, 4]:
            for csv_file in components_files:
                spreadsheet_path = os.path.join(
                    self.BASE_PATH, f"Channel {channel_id}", csv_file + ".csv"
                )
                pd.read_csv(spreadsheet_path)

    # # -----------------------Test SPARC4 operation mode ---------------------------------------

    def test_verify_sparc4_operation_mode_photometric(self):
        channel = Channel(self.CHANNEL)
        channel.write_sparc4_operation_mode("photometry")

    def test_verify_sparc4_operation_mode_error(self):
        channel = Channel(self.CHANNEL)
        with pytest.raises(ValueError):
            channel.write_sparc4_operation_mode("error")

    def test_verify_sparc4_operation_mode_polarimetric(self):
        channel = Channel(self.CHANNEL)
        channel.write_sparc4_operation_mode("polarimetry", "", "half")

    def test_verify_sparc4_operation_mode_polarimetric_polarizer(self):
        channel = Channel(self.CHANNEL)
        channel.write_sparc4_operation_mode("polarimetry", "polarizer", "half")

    def test_verify_sparc4_operation_mode_polarimetric_depolarizer(self):
        channel = Channel(self.CHANNEL)
        channel.write_sparc4_operation_mode("polarimetry", "depolarizer", "half")

    def test_verify_sparc4_operation_mode_polarimetric_quarter(self):
        channel = Channel(self.CHANNEL)
        channel.write_sparc4_operation_mode("polarimetry", "", "quarter")

    def test_verify_sparc4_operation_mode_polarimetric_wrong_retarder_value(self):
        channel = Channel(self.CHANNEL)
        with pytest.raises(ValueError):
            channel.write_sparc4_operation_mode("polarimetry", "", "")

    def test_verify_sparc4_operation_mode_polarimetric_wrong_calibration_wheel_value(
        self,
    ):
        channel = Channel(self.CHANNEL)
        with pytest.raises(ValueError):
            channel.write_sparc4_operation_mode("polarimetry", "error", "half")

    # # ---------------------------------miscelaneous----------------------------------------------

    def test_interpolate_spectral_response(self):
        channel = self.channel
        channel.obj_wavelength = self.OBJ_WAVELENGTH
        new_spectral_response = channel._interpolate_spectral_response(
            self.wavelength, self.spectral_response
        )
        assert np.allclose(self.interpolated_spectral_response, new_spectral_response)

    def test_get_spectral_response(self):
        new_spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "polarizer.csv"
        )
        assert np.allclose(
            new_spectral_response, self.interpolated_spectral_response, rtol=0.005
        )

    def test_get_spectral_response_custom(self):
        channel = self.channel
        csv_file = "polarizer_contrast_ratio.csv"
        channel.obj_wavelength = self.OBJ_WAVELENGTH
        new_contrast_ratio = channel._get_spectral_response_custom(
            csv_file, "Contrast ratio"
        )

        csv_file_name = os.path.join(self.BASE_PATH, csv_file)
        ss = pd.read_csv(csv_file_name)
        sys_wavelength, contrast_ratio = ss["Wavelength (nm)"], ss["Contrast ratio"]

        contrast_ratio = channel._interpolate_spectral_response(
            sys_wavelength, contrast_ratio
        )

        assert np.allclose(contrast_ratio, new_contrast_ratio)

    # ----------------------------------------------------------------------------------------------------

    def test_apply_ideal_polarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "polarizer.csv"
        )
        self.channel.write_sparc4_operation_mode("polarimetry", "ideal-polarizer")
        self.channel.sed = self.SED
        sed = np.multiply(spectral_response, self.SED)
        polarizer_matrix = calculate_polarizer_matrix(self.POLARIZER_ANGLE)
        sed = apply_matrix(polarizer_matrix, sed)

        self.channel._apply_ideal_polarizer(spectral_response)
        assert np.allclose(sed, self.channel.sed)

    def test_apply_real_polarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "polarizer.csv"
        )
        contrast_ratio = self.channel._get_spectral_response_custom(
            "polarizer_contrast_ratio.csv", "Contrast ratio"
        )
        sed = self.SED
        for idx, transmission in enumerate(spectral_response):
            contrast = contrast_ratio[idx]
            polarizer_matrix = calculate_polarizer_matrix(self.POLARIZER_ANGLE)
            sed[:, idx] = transmission * np.transpose(polarizer_matrix.dot(sed[:, idx]))
            sed[1, idx] *= 1 - 1 / contrast

        self.channel.write_sparc4_operation_mode("polarimetry", "polarizer")
        self.channel.sed = self.SED
        self.channel._apply_real_polarizer(spectral_response)
        assert np.allclose(sed, self.channel.sed)

    def test_apply_polarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "polarizer.csv"
        )
        self.channel.sed = self.SED
        self.channel._apply_real_polarizer(spectral_response)
        sed = self.channel.sed

        self.channel.sed = self.SED
        self.channel.obj_wavelength = self.OBJ_WAVELENGTH
        self.channel.write_sparc4_operation_mode("polarimetry", "polarizer")
        self.channel._apply_polarizer()

        assert np.allclose(sed, self.channel.sed)

    def test_apply_ideal_depolarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "depolarizer.csv"
        )

        tmp = self.SED
        sed = np.zeros((4, self.SED.shape[1]))
        sed[0] = tmp[0] * spectral_response

        self.channel.sed = self.SED
        self.channel._apply_ideal_depolarizer(spectral_response)
        assert np.allclose(sed, self.channel.sed)

    def test_apply_real_depolarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "depolarizer.csv"
        )
        sed = np.multiply(spectral_response, self.SED)
        for idx, wavelength in enumerate(self.OBJ_WAVELENGTH):
            depolarizer_matrix = calculate_depolarizer_matrix(wavelength)
            sed[:, idx] = np.transpose(depolarizer_matrix.dot(sed[:, idx]))

        self.channel.sed = self.SED
        self.channel._apply_real_depolarizer(spectral_response)
        assert np.allclose(sed, self.channel.sed)

    def test_apply_depolarizer(self):
        spectral_response = self.channel.get_spectral_response(
            self.OBJ_WAVELENGTH, "depolarizer.csv"
        )
        self.channel.sed = self.SED
        self.channel._apply_real_depolarizer(spectral_response)
        sed = self.channel.sed

        self.channel.sed = self.SED
        self.channel.obj_wavelength = self.OBJ_WAVELENGTH
        self.channel.write_sparc4_operation_mode("polarimetry", "depolarizer")
        self.channel._apply_depolarizer()

        assert np.allclose(sed, self.channel.sed)

    def test_apply_ideal_retarder_waveplate_half(self):
        RETARDER_MATRIX = calculate_retarder_matrix(180, 0)
        sed = apply_matrix(RETARDER_MATRIX, self.SED)

        self.channel.sed = self.SED
        self.channel.write_sparc4_operation_mode("polarimetry")
        self.channel._apply_ideal_waveplate()

        assert np.allclose(sed, self.channel.sed)

    def test_apply_real_retarder_waveplate_half(self):
        self.channel.obj_wavelength = self.OBJ_WAVELENGTH
        phase_difference = (
            self.channel._get_spectral_response_custom(
                "retarder_phase_diff_half.csv", "Retardance"
            )
            * 360
        )
        sed = self.SED
        for idx, phase in enumerate(phase_difference):
            RETARDER_MATRIX = calculate_retarder_matrix(phase, 0)
            sed[:, idx] = np.transpose(RETARDER_MATRIX.dot(sed[:, idx]))

        self.channel.write_sparc4_operation_mode("polarimetry")
        self.channel.sed = self.SED
        self.channel._apply_real_waveplate()

        assert np.allclose(sed, self.channel.sed)

    def test_apply_retarder_waveplate(self):
        self.channel.sed = self.SED
        self.channel.obj_wavelength = self.OBJ_WAVELENGTH
        self.channel.write_sparc4_operation_mode("polarimetry")
        self.channel._apply_real_waveplate()
        sed = self.channel.sed

        self.channel.sed = self.SED
        self.channel._apply_retarder_waveplate()

        assert np.allclose(sed, self.channel.sed)

    def test_apply_analyzer(self):
        ORD_RAY_MATRIX = calculate_polarizer_matrix(0)
        temp_1 = apply_matrix(ORD_RAY_MATRIX, self.SED)
        EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(90)
        temp_2 = apply_matrix(EXTRA_ORD_RAY_MATRIX, self.SED)

        self.channel.sed = self.SED
        self.channel._apply_analyzer()

        assert np.allclose(self.channel.sed, np.stack((temp_1[0], temp_2[0])))

    # ----------------------------------------------------------------------------------------------------

    def test_apply_phot_spectral_response(self):
        sed = self.SED
        for csv_file in self.PHOT_OPTICAL_COMPONENTS.values():
            if csv_file != "collimator.csv":
                csv_file = os.path.join(f"Channel {self.CHANNEL}", csv_file)
            spectral_response = self.channel.get_spectral_response(
                self.OBJ_WAVELENGTH, csv_file
            )
            sed = np.multiply(spectral_response, sed)

        self.channel.sed = self.SED
        self.channel._apply_photometric_spectral_response()

        assert np.allclose(self.channel.sed, sed)

    def test_apply_pol_spectral_response(self):
        sed = self.SED
        for csv_file in self.POL_OPTICAL_COMPONENTS.values():
            spectral_response = self.channel.get_spectral_response(
                self.OBJ_WAVELENGTH, csv_file
            )
            sed = np.multiply(spectral_response, sed)

        self.channel.sed = sed
        self.channel.write_sparc4_operation_mode("polarimetry", "polarizer")
        self.channel._apply_calibration_wheel()
        self.channel._apply_retarder_waveplate()
        self.channel._apply_analyzer()
        sed = self.channel.sed

        self.channel.sed = self.SED
        self.channel.obj_wavelength = self.OBJ_WAVELENGTH
        self.channel._apply_polarimetric_spectral_response()

        assert np.allclose(sed, self.channel.sed)
