# # -*- coding: utf-8 -*-
# Tests of the Channel Class
# Oct 17th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Response import Channel
from AIS.Spectral_Response._utils import calculate_polarizer_matrix, calculate_retarder_matrix
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep, interp1d
from math import atan, sqrt


obj_wavelength = np.linspace(400, 1100, 100)
CHANNEL = 1
BASE_PATH = os.path.join("AIS", "Spectral_Response", "channel")
_POL_OPTICAL_COMPONENTS = {"retarder": "retarder.csv",
                           "analyzer": "analyzer.csv",
                           }

_PHOT_OPTICAL_COMPONENTS = {
    "collimator": "collimator.csv",
    "dichroic": "dichroic.csv",
    "camera": "camera.csv",
    "ccd": "ccd.csv",
}


@pytest.fixture
def channel():
    return Channel(CHANNEL)


def _calc_polarizer_matrix(pol_angle: float):
    POLARIZER_MATRIX = 0.5 * np.asarray(
        [
            [1, cos(2 * pol_angle), sin(2 * pol_angle), 0],
            [
                cos(2 * pol_angle),
                cos(2 * pol_angle) ** 2,
                cos(2 * pol_angle) * sin(2 * pol_angle),
                0,
            ],
            [
                sin(2 * pol_angle),
                cos(2 * pol_angle) * sin(2 * pol_angle),
                sin(2 * pol_angle) ** 2,
                0,
            ],
            [0, 0, 0, 0],
        ]
    )
    return POLARIZER_MATRIX


def _calc_retarder_matrix(phase_difference, ret_angle):
    phase_difference = np.radians(phase_difference)
    ret_angle = np.radians(ret_angle)
    RETARDER_MATRIX = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * ret_angle) ** 2
                + sin(2 * ret_angle) ** 2 * cos(phase_difference),
                cos(2 * ret_angle)
                * sin(2 * ret_angle)
                * (1 - cos(phase_difference)),
                -sin(2 * ret_angle) * sin(phase_difference),
            ],
            [
                0,
                cos(2 * ret_angle)
                + sin(2 * ret_angle) * (1 - cos(phase_difference)),
                sin(2 * ret_angle) ** 2
                + cos(2 * ret_angle) ** 2 * cos(phase_difference),
                cos(2 * ret_angle) * sin(phase_difference),
            ],
            [
                0,
                sin(2 * ret_angle) * sin(phase_difference),
                -cos(2 * ret_angle) * sin(phase_difference),
                cos(phase_difference),
            ],
        ]
    )

    return RETARDER_MATRIX
# ---------------------------------------------------------------------------------------------------------------


def test_csv_file_name_pol(channel):
    assert channel._POL_OPTICAL_COMPONENTS == _POL_OPTICAL_COMPONENTS


def test_csv_file_name_phot(channel):
    assert channel._PHOT_OPTICAL_COMPONENTS == _PHOT_OPTICAL_COMPONENTS


def test_base_path(channel):
    assert channel._BASE_PATH == BASE_PATH


# --------------------------------------------------------------------------------------------

polarizer_path = os.path.join(BASE_PATH, "polarizer.csv")
ss = pd.read_csv(polarizer_path)
wv_polarizer = ss["Wavelength (nm)"]
sr_polarizer = ss["Transmitance (%)"] / 100


def test_read_csv_file_polarizer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(polarizer_path)
    assert np.allclose(new_wavelength, wv_polarizer)
    assert np.allclose(new_transmitance, sr_polarizer)


depolarizer_path = os.path.join(
    BASE_PATH, "depolarizer.csv")
ss = pd.read_csv(depolarizer_path)
wv_depolarizer = ss["Wavelength (nm)"]
sr_depolarizer = ss["Transmitance (%)"] / 100


def test_read_csv_file_depolarizer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(depolarizer_path)
    assert np.allclose(new_wavelength, wv_depolarizer)
    assert np.allclose(new_transmitance, sr_depolarizer)


retarder_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["retarder"])
ss = pd.read_csv(retarder_path)
wv_retarder = ss["Wavelength (nm)"]
sr_retarder = ss["Transmitance (%)"] / 100


def test_read_csv_file_retarder(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(retarder_path)
    assert np.allclose(new_wavelength, wv_retarder)
    assert np.allclose(new_transmitance, sr_retarder)


analyzer_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["analyzer"])
ss = pd.read_csv(analyzer_path)
wv_analyzer = ss["Wavelength (nm)"]
sr_analyzer = ss["Transmitance (%)"] / 100


def test_read_csv_file_analyzer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(analyzer_path)
    assert np.allclose(new_wavelength, wv_analyzer)
    assert np.allclose(new_transmitance, sr_analyzer)


collimator_path = os.path.join(
    BASE_PATH, _PHOT_OPTICAL_COMPONENTS["collimator"])
ss = pd.read_csv(collimator_path)
wv_collimator = ss["Wavelength (nm)"]
sr_collimator = ss["Transmitance (%)"] / 100


def test_read_csv_file_collimator(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(collimator_path)
    assert np.allclose(new_wavelength, wv_collimator)
    assert np.allclose(new_transmitance, sr_collimator)


dichroic_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["dichroic"]
)
ss = pd.read_csv(dichroic_path)
wv_dichroic = ss["Wavelength (nm)"]
sr_dichroic = ss["Transmitance (%)"] / 100


def test_read_csv_file_dichroic(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(dichroic_path)
    assert np.allclose(new_wavelength, wv_dichroic)
    assert np.allclose(new_transmitance, sr_dichroic)


camera_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["camera"])
ss = pd.read_csv(camera_path)
wv_camera = ss["Wavelength (nm)"]
sr_camera = ss["Transmitance (%)"] / 100


def test_read_csv_file_camera(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(camera_path)
    assert np.allclose(new_wavelength, wv_camera)
    assert np.allclose(new_transmitance, sr_camera)


ccd_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["ccd"])
ss = pd.read_csv(ccd_path)
wv_ccd = ss["Wavelength (nm)"]
sr_ccd = ss["Transmitance (%)"] / 100


def test_read_csv_file_ccd(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(ccd_path)
    assert np.allclose(new_wavelength, wv_ccd)
    assert np.allclose(new_transmitance, sr_ccd)

# -----------------------Test SPARC4 operation mode ---------------------------------------


def test_verify_sparc4_operation_mode_photometric():
    channel = Channel(CHANNEL)
    channel.write_sparc4_operation_mode("photometry")


def test_verify_sparc4_operation_mode_error():
    channel = Channel(CHANNEL)
    with pytest.raises(ValueError):
        channel.write_sparc4_operation_mode('error')


def test_verify_sparc4_operation_mode_polarimetric():
    channel = Channel(CHANNEL)
    channel.write_sparc4_operation_mode("polarimetry", '', 'half')


def test_verify_sparc4_operation_mode_polarimetric_polarizer():
    channel = Channel(CHANNEL)
    channel.write_sparc4_operation_mode("polarimetry", 'polarizer', 'half')


def test_verify_sparc4_operation_mode_polarimetric_depolarizer():
    channel = Channel(CHANNEL)
    channel.write_sparc4_operation_mode("polarimetry", 'depolarizer', 'half')


def test_verify_sparc4_operation_mode_polarimetric_quarter():
    channel = Channel(CHANNEL)
    channel.write_sparc4_operation_mode("polarimetry", '', 'quarter')


def test_verify_sparc4_operation_mode_polarimetric_wrong_retarder_value():
    channel = Channel(CHANNEL)
    with pytest.raises(ValueError):
        channel.write_sparc4_operation_mode("polarimetry", '', '')


def test_verify_sparc4_operation_mode_polarimetric_wrong_calibration_wheel_value():
    channel = Channel(CHANNEL)
    with pytest.raises(ValueError):
        channel.write_sparc4_operation_mode("polarimetry", 'error', 'half')


# ---------------------------------miscelaneous----------------------------------------------

# spl = splrep(wv_collimator, sr_collimator)
# new_spectral_response = splev(obj_wavelength, spl)

spl = interp1d(wv_collimator, sr_collimator,
               bounds_error=False, fill_value=0, kind='cubic')
new_spectral_response = spl(obj_wavelength)


def test_interpolate_spectral_response(channel):
    class_spectral_response = channel._interpolate_spectral_response(
        wv_collimator, sr_collimator, obj_wavelength
    )
    assert np.allclose(class_spectral_response, new_spectral_response)


def test_get_spectral_response(channel):
    class_spectral_response = channel.get_spectral_response(
        obj_wavelength, _PHOT_OPTICAL_COMPONENTS["collimator"]
    )
    assert np.allclose(class_spectral_response,
                       new_spectral_response, rtol=0.005)


def test_get_polarizer_contrast_ratio(channel):
    contrast_ratio = channel._get_polarizer_contrast_ratio(obj_wavelength)
    csv_file_name = os.path.join(
        'AIS', 'Spectral_Response', 'channel', 'polarizer_contrast_ratio.csv')

    ss = pd.read_csv(csv_file_name)
    sys_wavelength, new_contrast_ratio = np.asarray(
        ss["Wavelength (nm)"]), np.asarray(ss["Contrast ratio"])
    spl = interp1d(sys_wavelength, new_contrast_ratio,
                   bounds_error=False, fill_value=0, kind='cubic')
    assert np.allclose(contrast_ratio, spl(obj_wavelength))

    # ----------------------------------------------------------------------------------------------------


def test_resize_sed(channel):
    sed = np.ones(10)
    channel.sed = sed
    channel._resize_sed()
    temp = sed
    sed = np.zeros((4, 10))
    sed[0, :] = temp
    assert np.allclose(channel.sed, sed)


def test_resize_sed_2(channel):
    sed = np.ones((4, 10))
    channel.sed = sed
    channel._resize_sed()
    assert np.allclose(channel.sed, sed)


def test_calc_polarizer_matrix():
    pol_angle = 30
    radians = np.deg2rad(pol_angle)
    POLARIZER_MATRIX = 0.5 * np.asarray(
        [
            [1, cos(2 * radians), sin(2 * radians), 0],
            [
                cos(2 * radians),
                cos(2 * radians) ** 2,
                cos(2 * radians) * sin(2 * radians),
                0,
            ],
            [
                sin(2 * radians),
                cos(2 * radians) * sin(2 * radians),
                sin(2 * radians) ** 2,
                0,
            ],
            [0, 0, 0, 0],
        ]
    )
    polarizer_matrix = calculate_polarizer_matrix(pol_angle)
    assert np.allclose(polarizer_matrix, POLARIZER_MATRIX)


def test_calc_retarder_matrix():
    phase_diff = np.radians(30)
    pol_angle = np.radians(60)
    RETARDER_MATRIX = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * pol_angle) ** 2
                + sin(2 * pol_angle) ** 2 * cos(phase_diff),
                cos(2 * pol_angle)
                * sin(2 * pol_angle)
                * (1 - cos(phase_diff)),
                -sin(2 * pol_angle) * sin(phase_diff),
            ],
            [
                0,
                cos(2 * pol_angle)
                + sin(2 * pol_angle) * (1 - cos(phase_diff)),
                sin(2 * pol_angle) ** 2
                + cos(2 * pol_angle) ** 2 * cos(phase_diff),
                cos(2 * pol_angle) * sin(phase_diff),
            ],
            [
                0,
                sin(2 * pol_angle) * sin(phase_diff),
                -cos(2 * pol_angle) * sin(phase_diff),
                cos(phase_diff),
            ],
        ]
    )
    retarder_matrix = calculate_retarder_matrix(30, 60)
    assert np.allclose(retarder_matrix, RETARDER_MATRIX)


def test_apply_calibration_wheel_polarizer(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    spl = interp1d(wv_polarizer, sr_polarizer,
                   bounds_error=False, fill_value=0, kind='cubic')
    reduced_sed = np.multiply(spl(obj_wavelength), sed)

    contrast_ratio = channel._get_polarizer_contrast_ratio(obj_wavelength)
    for idx, value in enumerate(contrast_ratio):
        theta = np.rad2deg(atan(sqrt(value)))
        polarizer_matrix = calculate_polarizer_matrix(theta)
        reduced_sed[:, idx] = np.transpose(
            polarizer_matrix.dot(reduced_sed[:, idx]))

    channel.calibration_wheel = 'polarizer'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_calibration_wheel()
    assert np.allclose(channel.sed, reduced_sed, atol=1e-3)


def test_apply_calibration_wheel_depolarizer(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    spl = interp1d(wv_depolarizer, sr_depolarizer,
                   bounds_error=False, fill_value=0, kind='cubic')
    reduced_sed[0, :] = np.multiply(spl(obj_wavelength), reduced_sed[0])
    temp = reduced_sed
    reduced_sed = np.zeros((4, n))
    reduced_sed[0, :] = temp[0, :]

    channel.calibration_wheel = 'depolarizer'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_calibration_wheel()

    assert np.allclose(channel.sed, reduced_sed)


def test_apply_retarder_waveplate(channel):
    phase_diff = 180
    retarder_waveplate_angle = 0
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    RETARDER_MATRIX = calculate_retarder_matrix(
        phase_diff, retarder_waveplate_angle)
    reduced_sed = np.transpose([RETARDER_MATRIX.dot(sed[:, i])
                                for i in range(sed.shape[1])])

    channel.retarder_waveplate_angle = retarder_waveplate_angle
    channel.retarder_waveplate = 'half'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_retarder_waveplate()
    assert np.allclose(channel.sed, reduced_sed)


def test_apply_analyzer(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))

    ORD_RAY_MATRIX = calculate_polarizer_matrix(0)
    temp_1 = np.transpose([ORD_RAY_MATRIX.dot(sed[:, i])
                          for i in range(sed.shape[1])])

    EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(90)
    temp_2 = np.transpose([EXTRA_ORD_RAY_MATRIX.dot(sed[:, i])
                           for i in range(sed.shape[1])])

    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_analyzer()
    assert np.allclose(channel.sed, np.stack((temp_1[0], temp_2[0])))

# ----------------------------------------------------------------------------------------------------


def test_apply_phot_spectral_response(channel):
    n = len(obj_wavelength)
    sed = np.linspace(400, 1100, n)
    spl = interp1d(wv_collimator, sr_collimator,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_esd = np.multiply(sed, new_spectral_response)
    spl = interp1d(wv_dichroic, sr_dichroic,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = interp1d(wv_camera, sr_camera,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = interp1d(wv_ccd, sr_ccd,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_photometric_spectral_response()
    assert np.allclose(channel.sed, reduced_esd)


def test_apply_pol_spectral_response_1(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    phase_diff = 180
    retarder_waveplate_angle = 0

    spl = interp1d(wv_polarizer, sr_polarizer,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed = np.multiply(sed, new_spectral_response)
    csv_file_name = os.path.join(
        BASE_PATH, 'polarizer_contrast_ratio.csv')
    ss = pd.read_csv(csv_file_name)
    sys_wavelength, contrast_ratio = np.asarray(
        ss["Wavelength (nm)"]), np.asarray(ss["Contrast ratio"])
    spl = interp1d(sys_wavelength, contrast_ratio,
                   bounds_error=False, fill_value=0, kind='cubic')
    contrast_ratio = spl(obj_wavelength)
    for idx, value in enumerate(contrast_ratio):
        theta = np.rad2deg(atan(sqrt(value)))
        POLARIZER_MATRIX = calculate_polarizer_matrix(theta)
        reduced_sed[:, idx] = np.transpose(
            POLARIZER_MATRIX.dot(reduced_sed[:, idx]))

    spl = interp1d(wv_retarder, sr_retarder,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed = np.multiply(reduced_sed, new_spectral_response)
    RETARDER_MATRIX = calculate_retarder_matrix(
        phase_diff, retarder_waveplate_angle)
    reduced_sed = np.transpose([RETARDER_MATRIX.dot(reduced_sed[:, i])
                                for i in range(reduced_sed.shape[1])])

    spl = interp1d(wv_analyzer, sr_analyzer,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed = np.multiply(reduced_sed, new_spectral_response)
    ORD_RAY_MATRIX = calculate_polarizer_matrix(0)
    temp_1 = np.transpose([ORD_RAY_MATRIX.dot(reduced_sed[:, i])
                           for i in range(reduced_sed.shape[1])])
    EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(90)
    temp_2 = np.transpose([EXTRA_ORD_RAY_MATRIX.dot(reduced_sed[:, i])
                           for i in range(reduced_sed.shape[1])])

    channel.calibration_wheel = 'polarizer'
    channel.retarder_waveplate_angle = retarder_waveplate_angle
    channel.retarder_waveplate = 'half'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_polarimetric_spectral_response()
    assert np.allclose(channel.sed,
                       np.stack((temp_1[0], temp_2[0])))
