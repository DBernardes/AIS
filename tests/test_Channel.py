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
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep, interp1d


obj_wavelength = np.linspace(400, 1100, 100)
CHANNEL = 1
BASE_PATH = os.path.join("AIS", "Spectral_Response", "channel")
_POL_OPTICAL_COMPONENTS = {
    "polarizer": "polarizer.csv",
    "depolarizer": "depolarizer.csv",
    "retarder": "retarder.csv",
    "analyser": "analyser.csv",
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
    pol_angle = np.radians(pol_angle)
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

polarizer_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["polarizer"])
ss = pd.read_csv(polarizer_path)
wv_polarizer = ss["Wavelength (nm)"]
sr_polarizer = ss["Transmitance (%)"] / 100


def test_read_csv_file_polarizer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(polarizer_path)
    assert np.allclose(new_wavelength, wv_polarizer)
    assert np.allclose(new_transmitance, sr_polarizer)


depolarizer_path = os.path.join(
    BASE_PATH, _POL_OPTICAL_COMPONENTS["depolarizer"])
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


analyser_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["analyser"])
ss = pd.read_csv(analyser_path)
wv_analyser = ss["Wavelength (nm)"]
sr_analyser = ss["Transmitance (%)"] / 100


def test_read_csv_file_analyser(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(analyser_path)
    assert np.allclose(new_wavelength, wv_analyser)
    assert np.allclose(new_transmitance, sr_analyser)


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
        channel.write_sparc4_operation_mode("polarimetry", '', 'half')


# --------------------------------------------------------------------------------------------

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


def test_calc_polarizer_matrix(channel):
    pol_angle = np.rad2deg(pi/3)
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
    polarizer_matrix = channel._calc_polarizer_matrix(pol_angle)
    assert np.allclose(polarizer_matrix, POLARIZER_MATRIX)


def test_calc_retarder_matrix(channel):
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
    retarder_matrix = channel._calc_retarder_matrix(30, 60)
    assert np.allclose(retarder_matrix, RETARDER_MATRIX)


def test_apply_calibration_wheel_polarizer(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    spl = interp1d(wv_polarizer, sr_polarizer,
                   bounds_error=False, fill_value=0, kind='cubic')
    reduced_sed[0, :] = np.multiply(spl(obj_wavelength), sed[0])
    POLARIZER_MATRIX = _calc_polarizer_matrix(0)
    reduced_sed = np.transpose([POLARIZER_MATRIX.dot(reduced_sed[:, i])
                                for i in range(reduced_sed.shape[1])])
    channel.calibration_wheel = 'polarizer'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_calibration_wheel()
    assert np.allclose(channel.sed, reduced_sed)


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


# def test_apply_calibration_wheel_error(channel):
#     channel.calibration_wheel = 'error'
#     channel.obj_wavelength = obj_wavelength
#     with pytest.raises(ValueError):
#         channel._apply_calibration_wheel()


def test_apply_retarder_waveplate(channel):
    phase_diff = 90
    retarder_waveplate_angle = 0
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    spl = interp1d(wv_retarder, sr_retarder,
                   bounds_error=False, fill_value=0, kind='cubic')
    reduced_sed[0, :] = np.multiply(spl(obj_wavelength), reduced_sed[0])
    RETARDER_MATRIX = _calc_retarder_matrix(
        phase_diff, retarder_waveplate_angle)
    reduced_sed = np.transpose([RETARDER_MATRIX.dot(reduced_sed[:, i])
                                for i in range(reduced_sed.shape[1])])
    channel.retarder_waveplate_angle = retarder_waveplate_angle
    channel.retarder_waveplate = 'half'
    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_retarder_waveplate()
    assert np.allclose(channel.sed, reduced_sed)


def test_apply_analyser(channel):
    n = len(obj_wavelength)
    sed = np.ones((4, n))
    reduced_sed = np.ones((4, n))
    spl = interp1d(wv_analyser, sr_analyser,
                   bounds_error=False, fill_value=0, kind='cubic')
    reduced_sed[0, :] = np.multiply(spl(obj_wavelength), reduced_sed[0])

    ORD_RAY_MATRIX = _calc_polarizer_matrix(0)
    temp_1 = np.transpose([ORD_RAY_MATRIX.dot(reduced_sed[:, i])
                          for i in range(reduced_sed.shape[1])])

    EXTRA_ORD_RAY_MATRIX = _calc_polarizer_matrix(90)
    temp_2 = np.transpose([EXTRA_ORD_RAY_MATRIX.dot(reduced_sed[:, i])
                           for i in range(reduced_sed.shape[1])])

    channel.sed = sed
    channel.obj_wavelength = obj_wavelength
    channel._apply_analyser()
    assert np.allclose(channel.sed, np.stack((temp_1[0], temp_2[0])))

# ----------------------------------------------------------------------------------------------------


def test_apply_phot_spectral_response(channel):
    n = len(obj_wavelength)
    sed = np.linspace(100, 360, n)
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
    phase_diff = 90
    retarder_waveplate_angle = 0

    spl = interp1d(wv_polarizer, sr_polarizer,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed[0, :] = np.multiply(sed[0], new_spectral_response)
    POLARIZER_MATRIX = _calc_polarizer_matrix(0)
    reduced_sed = np.transpose([POLARIZER_MATRIX.dot(reduced_sed[:, i])
                                for i in range(reduced_sed.shape[1])])

    spl = interp1d(wv_retarder, sr_retarder,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed[0, :] = np.multiply(reduced_sed[0], new_spectral_response)
    RETARDER_MATRIX = _calc_retarder_matrix(
        phase_diff, retarder_waveplate_angle)
    reduced_sed = np.transpose([RETARDER_MATRIX.dot(reduced_sed[:, i])
                                for i in range(reduced_sed.shape[1])])

    spl = interp1d(wv_analyser, sr_analyser,
                   bounds_error=False, fill_value=0, kind='cubic')
    new_spectral_response = spl(obj_wavelength)
    reduced_sed[0, :] = np.multiply(reduced_sed[0], new_spectral_response)
    ORD_RAY_MATRIX = _calc_polarizer_matrix(0)
    temp_1 = np.transpose([ORD_RAY_MATRIX.dot(reduced_sed[:, i])
                           for i in range(reduced_sed.shape[1])])
    EXTRA_ORD_RAY_MATRIX = _calc_polarizer_matrix(90)
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


# def test_apply_spectral_response(channel):
#     n = len(obj_wavelength)
#     esd = np.linspace(100, 360, n)
#     spl = interp1d(wv_retarder, sr_retarder,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(esd, new_spectral_response)
#     spl = interp1d(wv_analyser, sr_analyser,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(reduced_esd, new_spectral_response)
#     spl = interp1d(wv_collimator, sr_collimator,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(reduced_esd, new_spectral_response)
#     spl = interp1d(wv_dichroic, sr_dichroic,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(reduced_esd, new_spectral_response)
#     spl = interp1d(wv_camera, sr_camera,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(reduced_esd, new_spectral_response)
#     spl = interp1d(wv_ccd, sr_ccd,
#                    bounds_error=False, fill_value=0, kind='cubic')
#     new_spectral_response = spl(obj_wavelength)
#     reduced_esd = np.multiply(reduced_esd, new_spectral_response)
#     channel.calibration_wheel = ''
#     channel.acquisition_mode = 'polarimetry'
#     class_reduced_esd = channel.apply_spectral_response(
#         esd, obj_wavelength)
#     assert np.allclose(class_reduced_esd, reduced_esd)
