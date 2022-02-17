# -*- coding: utf-8 -*-

"""
SPARC4 spectral response tests.

This script tests the operation of the SPARC4 spectral response classes.
"""

import os

import numpy as np
import pandas as pd
import pytest
from AIS.SPARC4_Spectral_Response import (
    Abstract_SPARC4_Spectral_Response,
    Concrete_SPARC4_Spectral_Response_1,
    Concrete_SPARC4_Spectral_Response_2,
    Concrete_SPARC4_Spectral_Response_3,
    Concrete_SPARC4_Spectral_Response_4,
)
from scipy.interpolate import splev, splrep

from .AIS_spectral_response_curves import (
    POLARIZER_90_MATRIX,
    POLARIZER_MATRIX,
    THETA_POL,
    analyser_transmitance,
    calculate_retarder_matrix,
    camera_c0,
    camera_c1,
    camera_c2,
    camera_c3,
    camera_c4,
    ccd_transmitance_c0,
    ccd_transmitance_c1,
    ccd_transmitance_c2,
    ccd_transmitance_c3,
    ccd_transmitance_c4,
    colimator_transmitance,
    depolarizer_transmitance,
    dichroic_c0,
    dichroic_c1,
    dichroic_c2,
    dichroic_c3,
    dichroic_c4,
    polarizer_transmitance,
    retardance_half,
    retardance_quarter,
    retarder_transmitance,
    star_specific_photons_per_second,
    wavelength_interval,
    wavelength_interval_len,
)

# -------------------------------------------------------------------------------------------------------------


@pytest.fixture
def abs_s4_sr():
    chc = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    return chc


@pytest.fixture
def c1_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_1(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    return chc


@pytest.fixture
def c2_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_2(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    return chc


@pytest.fixture
def c3_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_3(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    return chc


@pytest.fixture
def c4_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_4(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    return chc


def multiply_matrices(matrix, specific_photons_per_second):
    new_specific_photons_per_second = []
    for array in specific_photons_per_second:
        for idx, value in enumerate(array[0]):
            array[:, idx] = np.dot(matrix, array[:, idx])
        new_specific_photons_per_second.append(array)
    return new_specific_photons_per_second


def apply_optical_component_transmission(
    component_transmission, specific_photons_per_second
):
    new_specific_photons_per_second = []
    for array in specific_photons_per_second:
        array[0] = np.multiply(array[0], component_transmission)
        new_specific_photons_per_second.append(array)
    return new_specific_photons_per_second


# -------------------- Initialize the class -----------------------


def test_specific_photons_per_second(abs_s4_sr):
    specific_photons_per_second = abs_s4_sr.read_specific_photons_per_second()
    assert np.allclose(specific_photons_per_second, star_specific_photons_per_second)


def test_specific_photons_per_second_c1(c1_s4_sr):
    c1_specific_photons_per_second = c1_s4_sr.read_specific_photons_per_second()
    assert np.allclose(c1_specific_photons_per_second, star_specific_photons_per_second)


# -------------------- Channel ID -----------------------


def test_channel_ID_abs(abs_s4_sr):
    assert abs_s4_sr.get_channel_ID() == 0


def test_channel_ID_c1(c1_s4_sr):
    assert c1_s4_sr.get_channel_ID() == 1


def test_channel_ID_c2(c2_s4_sr):
    assert c2_s4_sr.get_channel_ID() == 2


def test_channel_ID_c3(c3_s4_sr):
    assert c3_s4_sr.get_channel_ID() == 3


def test_channel_ID_c4(c4_s4_sr):
    assert c4_s4_sr.get_channel_ID() == 4


# -------------------- Apply spectral response  -----------------------


def test_calibration_wheel_polarizer(abs_s4_sr):
    abs_s4_sr.apply_calibration_wheel("polarizer")
    new_star_specific_photons_per_second = multiply_matrices(
        POLARIZER_MATRIX, star_specific_photons_per_second
    )
    new_star_specific_photons_per_second = apply_optical_component_transmission(
        polarizer_transmitance, new_star_specific_photons_per_second
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_calibration_wheel_depolarizer(abs_s4_sr):
    abs_s4_sr.apply_calibration_wheel("depolarizer")
    new_star_specific_photons_per_second = [np.zeros((4, wavelength_interval_len))]
    new_star_specific_photons_per_second[0][0] = star_specific_photons_per_second[0][0]
    new_star_specific_photons_per_second = apply_optical_component_transmission(
        depolarizer_transmitance, new_star_specific_photons_per_second
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_calibration_wheel_empty(abs_s4_sr):
    abs_s4_sr.apply_calibration_wheel("empty")
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        star_specific_photons_per_second,
    )


def test_retarder_quarter(abs_s4_sr):
    abs_s4_sr.apply_retarder("quarter")
    new_star_specific_photons_per_second = star_specific_photons_per_second
    for idx, value in enumerate(retardance_quarter):
        retarder_matrix = calculate_retarder_matrix(value, THETA_POL)
        new_star_specific_photons_per_second[0][:, idx] = np.dot(
            retarder_matrix, new_star_specific_photons_per_second[0][:, idx]
        )
    new_star_specific_photons_per_second[0][0] = np.multiply(
        retarder_transmitance, new_star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_retarder_half(abs_s4_sr):
    abs_s4_sr.apply_retarder("half")
    new_star_specific_photons_per_second = star_specific_photons_per_second
    for idx, value in enumerate(retardance_half):
        retarder_matrix = calculate_retarder_matrix(value, THETA_POL)
        new_star_specific_photons_per_second[0][:, idx] = np.dot(
            retarder_matrix, new_star_specific_photons_per_second[0][:, idx]
        )
    new_star_specific_photons_per_second[0][0] = np.multiply(
        retarder_transmitance, new_star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_analyzer(abs_s4_sr):
    abs_s4_sr.apply_analyser()

    ordinary_ray = multiply_matrices(POLARIZER_MATRIX, star_specific_photons_per_second)
    ordinary_ray = apply_optical_component_transmission(
        analyser_transmitance, ordinary_ray
    )[0]

    extra_ordinary_ray = multiply_matrices(
        POLARIZER_90_MATRIX, star_specific_photons_per_second
    )
    extra_ordinary_ray = apply_optical_component_transmission(
        analyser_transmitance, extra_ordinary_ray
    )[0]
    new_star_specific_photons_per_second = [ordinary_ray, extra_ordinary_ray]
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_collimator(abs_s4_sr):
    abs_s4_sr.apply_collimator()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        retarder_transmitance, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_dichroic_abs(abs_s4_sr):
    abs_s4_sr.apply_dichroic()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        dichroic_c0, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_dichroic_c1(c1_s4_sr):
    c1_s4_sr.apply_dichroic()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        dichroic_c1, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c1_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_dichroic_c2(c2_s4_sr):
    c2_s4_sr.apply_dichroic()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        dichroic_c2, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c2_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_dichroic_c3(c3_s4_sr):
    c3_s4_sr.apply_dichroic()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        dichroic_c3, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c3_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_dichroic_c4(c4_s4_sr):
    c4_s4_sr.apply_dichroic()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        dichroic_c4, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c4_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_camera_abs(abs_s4_sr):
    abs_s4_sr.apply_camera()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        camera_c0, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_camera_c1(c1_s4_sr):
    c1_s4_sr.apply_camera()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        camera_c1, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c1_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_camera_c2(c2_s4_sr):
    c2_s4_sr.apply_camera()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        camera_c2, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c2_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_camera_c3(c3_s4_sr):
    c3_s4_sr.apply_camera()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        camera_c3, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c3_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_camera_c4(c4_s4_sr):
    c4_s4_sr.apply_camera()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        camera_c4, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c4_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_ccd_abs(abs_s4_sr):
    abs_s4_sr.apply_ccd()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        ccd_transmitance_c0, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_ccd_c1(c1_s4_sr):
    c1_s4_sr.apply_ccd()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        ccd_transmitance_c1, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c1_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_ccd_c2(c2_s4_sr):
    c2_s4_sr.apply_ccd()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        ccd_transmitance_c2, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c2_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_ccd_c3(c3_s4_sr):
    c3_s4_sr.apply_ccd()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        ccd_transmitance_c3, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c3_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


def test_ccd_c4(c4_s4_sr):
    c4_s4_sr.apply_ccd()
    new_star_specific_photons_per_second = star_specific_photons_per_second
    new_star_specific_photons_per_second[0][0] = np.multiply(
        ccd_transmitance_c4, star_specific_photons_per_second[0][0]
    )
    assert np.allclose(
        c4_s4_sr.read_specific_photons_per_second(),
        new_star_specific_photons_per_second,
    )


# --------------------write star_specific_photons_per_second--------------------


def test_write_specific_flux():
    star_specific_photons_per_second = np.asarray(
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
    )
    wavelength_interval = range(400, 1150, 50)
    s4_sr = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    s4_sr.write_specific_photons_per_second(star_specific_photons_per_second)
    assert np.allclose(
        s4_sr.read_specific_photons_per_second(), star_specific_photons_per_second
    )


# ----------------------- read_spreadsheet---------------------------

path = os.path.join("AIS", "SPARC4_Spectral_Response")


def test_read_spreadsheet_collimator(abs_s4_sr):
    file = os.path.join(path, "collimator.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_1(abs_s4_sr):
    file = os.path.join(path, "Channel 0", "dichroic_1.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_2(abs_s4_sr):
    file = os.path.join(path, "Channel 0", "dichroic_2.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera(abs_s4_sr):
    file = os.path.join(path, "Channel 0", "camera.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd(abs_s4_sr):
    file = os.path.join(path, "Channel 0", "ccd.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_1_1(abs_s4_sr):
    file = os.path.join(path, "Channel 1", "dichroic_1.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_1_2(abs_s4_sr):
    file = os.path.join(path, "Channel 1", "dichroic_2.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_1(abs_s4_sr):
    file = os.path.join(path, "Channel 1", "camera.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_1(abs_s4_sr):
    file = os.path.join(path, "Channel 1", "ccd.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_2_1(abs_s4_sr):
    file = os.path.join(path, "Channel 2", "dichroic_1.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_2_2(abs_s4_sr):
    file = os.path.join(path, "Channel 2", "dichroic_2.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_2(abs_s4_sr):
    file = os.path.join(path, "Channel 2", "camera.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_2(abs_s4_sr):
    file = os.path.join(path, "Channel 2", "ccd.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_3_1(abs_s4_sr):
    file = os.path.join(path, "Channel 3", "dichroic_1.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_3_2(abs_s4_sr):
    file = os.path.join(path, "Channel 3", "dichroic_2.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_3(abs_s4_sr):
    file = os.path.join(path, "Channel 3", "camera.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_3(abs_s4_sr):
    file = os.path.join(path, "Channel 3", "ccd.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_4_1(abs_s4_sr):
    file = os.path.join(path, "Channel 4", "dichroic_1.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_4_2(abs_s4_sr):
    file = os.path.join(path, "Channel 4", "dichroic_2.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_4(abs_s4_sr):
    file = os.path.join(path, "Channel 4", "camera.csv")
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_4(abs_s4_sr):
    file = os.path.join(path, "Channel 4", "ccd.csv")
    abs_s4_sr._read_spreadsheet(file)


# ----------------------- miscelaneous ----------------------------


def test_multiply_matrices(abs_s4_sr):
    a = np.ones((4, 4))
    abs_s4_sr._multiply_matrices(a)
    assert np.allclose(
        abs_s4_sr.read_specific_photons_per_second(), star_specific_photons_per_second
    )


def test_calculate_spline():
    transmitance = np.ones((1, wavelength_interval_len))[0]
    chc = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    chc.write_specific_photons_per_second(star_specific_photons_per_second)
    new_transmitance = chc._calculate_spline(transmitance, wavelength_interval)
    assert np.allclose(new_transmitance, transmitance)


def test_apply_optical_component_transmission(abs_s4_sr):
    file = os.path.join("AIS", "SPARC4_Spectral_Response", "polarizer.csv")
    abs_s4_sr._apply_optical_component_transmission(file)
    ss = pd.read_csv(file, skiprows=1, decimal=".", dtype=np.float64)
    component_wavelength = ss["(nm)"]
    component_transmitance = ss["(%)"]
    spl = splrep(component_wavelength, component_transmitance)
    component_transmitance = splev(wavelength_interval, spl)
    new_star_specific_photons_per_second = []
    for array in star_specific_photons_per_second:
        array[0] = np.multiply(array[0], component_transmitance)
        new_star_specific_photons_per_second.append(array)
    assert np.allclose(
        new_star_specific_photons_per_second,
        abs_s4_sr.read_specific_photons_per_second(),
    )


def test_calculate_retarder_matrix(abs_s4_sr):
    retarder_matrix = calculate_retarder_matrix(90, THETA_POL)
    new_retarder_matrix = abs_s4_sr._calculate_retarder_matrix(90)
    assert np.allclose(retarder_matrix, new_retarder_matrix)
