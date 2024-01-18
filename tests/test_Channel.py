# # -*- coding: utf-8 -*-
# Tests of the Channel Class
# Oct 17th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest, unittest
from AIS.Spectral_Response import Channel
from AIS.Spectral_Response._utils import (
    calculate_polarizer_matrix,
    calculate_retarder_matrix,
    calculate_depolarizer_matrix,
)
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep, interp1d, PchipInterpolator
from math import atan, sqrt
from copy import copy


class Test_Channel(unittest.TestCase):
    CHANNEL = 1
    OBJ_WAVELENGTH = np.linspace(400, 1100, 100)
    SED = np.ones(100)
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

    def test_base_path(self):
        assert self.channel.BASE_PATH == self.BASE_PATH

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

    #     # ----------------------------------------------------------------------------------------------------

    # def test_apply_calibration_wheel_polarizer(channel):
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     spl = interp1d(
    #         wv_polarizer, sr_polarizer, bounds_error=False, fill_value=0, kind="cubic"
    #     )

    #     channel.calibration_wheel = "polarizer"
    #     channel.sed = copy(sed)
    #     channel._apply_calibration_wheel()

    #     spl = splrep(
    #         wv_polarizer,
    #         sr_polarizer,
    #     )
    #     spectral_response = splev(obj_wavelength, spl)
    #     spectral_response[spectral_response < 0] = 0

    #     contrast_ratio = channel._get_spectral_response_custom(
    #         "polarizer_contrast_ratio.csv", "Contrast ratio"
    #     )
    #     for idx, transmission in enumerate(spectral_response):
    #         contrast = contrast_ratio[idx]
    #         polarizer_matrix = calculate_polarizer_matrix(POLARIZER_ANGLE)
    #         sed[:, idx] = transmission * np.transpose(polarizer_matrix.dot(sed[:, idx]))
    #         sed[1, idx] *= 1 - 1 / contrast

    #     assert np.allclose(channel.sed, sed, atol=1e-3)

    # def test_apply_calibration_wheel_ideal_polarizer(channel):
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     spl = interp1d(
    #         wv_polarizer, sr_polarizer, bounds_error=False, fill_value=0, kind="cubic"
    #     )

    #     channel.calibration_wheel = "ideal-polarizer"
    #     channel.sed = copy(sed)
    #     channel._apply_calibration_wheel()

    #     spectral_response = channel.get_spectral_response(obj_wavelength, "polarizer.csv")
    #     sed = np.multiply(spectral_response, sed)
    #     polarizer_matrix = calculate_polarizer_matrix(0)
    #     sed = _apply_matrix(polarizer_matrix, sed)

    #     assert np.allclose(sed, channel.sed)

    # def test_apply_calibration_wheel_ideal_depolarizer(channel):
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     channel.calibration_wheel = "ideal-depolarizer"
    #     channel.sed = copy(sed)
    #     channel._apply_calibration_wheel()

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_depolarizer, sr_depolarizer, obj_wavelength
    #     )
    #     sed = np.multiply(new_spectral_response, sed)
    #     reduced_sed = np.zeros((4, n))
    #     reduced_sed[0] = sed[0]

    #     assert np.allclose(channel.sed, reduced_sed, atol=1e-3)

    # def test_apply_calibration_wheel_depolarizer(channel):
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     channel.calibration_wheel = "depolarizer"
    #     channel.sed = copy(sed)
    #     channel._apply_calibration_wheel()

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_depolarizer, sr_depolarizer, obj_wavelength
    #     )
    #     sed = np.multiply(new_spectral_response, sed)
    #     for idx, wavelength in enumerate(obj_wavelength):
    #         depolarizer_matrix = calculate_depolarizer_matrix(wavelength)
    #         sed[:, idx] = np.transpose(depolarizer_matrix.dot(sed[:, idx]))

    #     assert np.allclose(channel.sed, sed, atol=1e-3)

    # def test_apply_retarder_waveplate(channel):
    #     phase_diff = "half"
    #     retarder_waveplate_angle = 0
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     channel.retarder_waveplate_angle = retarder_waveplate_angle
    #     channel.retarder_waveplate = phase_diff
    #     channel.sed = copy(sed)
    #     channel._apply_retarder_waveplate()

    #     phase_difference = (
    #         channel._get_spectral_response_custom(
    #             f"retarder_phase_diff_{phase_diff}.csv", "Retardance"
    #         )
    #         * 360
    #     )

    #     for idx, phase in enumerate(phase_difference):
    #         RETARDER_MATRIX = calculate_retarder_matrix(phase, retarder_waveplate_angle)
    #         sed[:, idx] = np.transpose(RETARDER_MATRIX.dot(sed[:, idx]))

    #     assert np.allclose(channel.sed, sed)

    # def test_apply_ideal_retarder_waveplate(channel):
    #     phase_diff = "ideal-half"
    #     retarder_waveplate_angle = 0
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     channel.retarder_waveplate_angle = retarder_waveplate_angle
    #     channel.retarder_waveplate = phase_diff
    #     channel.sed = copy(sed)
    #     channel._apply_retarder_waveplate()

    #     RETARDER_MATRIX = calculate_retarder_matrix(180, retarder_waveplate_angle)
    #     sed = _apply_matrix(RETARDER_MATRIX, sed)
    #     assert np.allclose(channel.sed, sed)

    # def test_apply_analyzer(channel):
    #     n = len(obj_wavelength)
    #     sed = np.ones((4, n))
    #     channel.sed = copy(sed)
    #     channel._apply_analyzer()

    #     ORD_RAY_MATRIX = calculate_polarizer_matrix(0)
    #     temp_1 = _apply_matrix(ORD_RAY_MATRIX, copy(sed))

    #     EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(90)
    #     temp_2 = _apply_matrix(EXTRA_ORD_RAY_MATRIX, copy(sed))

    #     assert np.allclose(channel.sed, np.stack((temp_1[0], temp_2[0])))

    # # ----------------------------------------------------------------------------------------------------

    # def test_apply_phot_spectral_response(channel):
    #     n = len(obj_wavelength)
    #     sed = np.zeros((4, n))
    #     sed = np.linspace(400, 1100, n)
    #     channel.sed = sed
    #     channel._apply_photometric_spectral_response()

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_collimator, sr_collimator, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)
    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_dichroic, sr_dichroic, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)
    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_camera, sr_camera, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)
    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_ccd, sr_ccd, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)

    #     assert np.allclose(channel.sed, sed)

    # def test_apply_pol_spectral_response(channel):
    #     n = len(obj_wavelength)
    #     sed = np.zeros((4, n))
    #     sed[0] = np.ones(n)
    #     phase_diff = "half"
    #     retarder_waveplate_angle = 0

    #     channel.calibration_wheel = "polarizer"
    #     channel.retarder_waveplate_angle = retarder_waveplate_angle
    #     channel.retarder_waveplate = phase_diff
    #     channel.sed = copy(sed)
    #     channel._apply_polarimetric_spectral_response()

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_polarizer, sr_polarizer, obj_wavelength
    #     )
    #     new_contrast_ratio = _interpolate_spectral_response(
    #         sys_wavelength, contrast_ratio, obj_wavelength
    #     )
    #     for idx, value in enumerate(new_contrast_ratio):
    #         theta = np.rad2deg(atan(1 / sqrt(value)))
    #         total_transmission = new_spectral_response[idx] * (1 + 1 / value)
    #         POLARIZER_MATRIX = calculate_polarizer_matrix(theta)
    #         sed[:, idx] = total_transmission * np.transpose(
    #             POLARIZER_MATRIX.dot(sed[:, idx])
    #         )

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_retarder, sr_retarder, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)
    #     phase_difference = channel._get_spectral_response_custom(
    #         f"retarder_phase_diff_{phase_diff}.csv", "Retardance"
    #     )
    #     for idx, phase in enumerate(phase_difference):
    #         phase = np.rad2deg(phase)
    #         RETARDER_MATRIX = calculate_retarder_matrix(phase, retarder_waveplate_angle)
    #         sed[:, idx] = np.transpose(RETARDER_MATRIX.dot(sed[:, idx]))

    #     new_spectral_response = _interpolate_spectral_response(
    #         wv_analyzer, sr_analyzer, obj_wavelength
    #     )
    #     sed = np.multiply(sed, new_spectral_response)
    #     ORD_RAY_MATRIX = calculate_polarizer_matrix(0)
    #     temp_1 = _apply_matrix(ORD_RAY_MATRIX, copy(sed))
    #     EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(90)
    #     temp_2 = _apply_matrix(EXTRA_ORD_RAY_MATRIX, copy(sed))

    #     assert np.allclose(channel.sed, np.stack((temp_1[0], temp_2[0])), atol=1e-3)
