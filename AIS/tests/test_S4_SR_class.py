# -*- coding: utf-8 -*-
"""SPARC4 spectral response tests.

This script tests the operation of the SPARC4 spectral response classes.
"""

import numpy as np
import pytest
from AIS.SPARC4_Spectral_Response import (
    Abstract_SPARC4_Spectral_Response,
    Concrete_SPARC4_Spectral_Response_1,
    Concrete_SPARC4_Spectral_Response_2,
    Concrete_SPARC4_Spectral_Response_3,
    Concrete_SPARC4_Spectral_Response_4,
)

wavelength_interval = range(350, 1100, 50)
n = len(wavelength_interval)
specific_flux = np.ones((4, n))


@pytest.fixture
def abs_s4_sr():
    chc = Abstract_SPARC4_Spectral_Response()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    return chc


@pytest.fixture
def c1_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_1()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    return chc


@pytest.fixture
def c2_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_2()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    return chc


@pytest.fixture
def c3_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_3()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    return chc


@pytest.fixture
def c4_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_4()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    return chc


# -------------------- Initialize the class -----------------------


def test_specific_flux_abs(abs_s4_sr):
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_specific_flux_c1(c1_s4_sr):
    vec = c1_s4_sr.get_specific_flux()[0, :]
    assert vec.all() == specific_flux.all()


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


def test_calibration_wheel(abs_s4_sr):
    abs_s4_sr.apply_polarimetric_component_spectral_response("calibration_wheel")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_retarder(abs_s4_sr):
    abs_s4_sr.apply_polarimetric_component_spectral_response("retarder")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_analyzer(abs_s4_sr):
    abs_s4_sr.apply_polarimetric_component_spectral_response("analyser")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_collimator(abs_s4_sr):
    abs_s4_sr.collimator()
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux[0, :].all()


def test_dichroic_abs(abs_s4_sr):
    abs_s4_sr.apply_photometric_component_spectral_response("dichroic")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_dichroic_c1(c1_s4_sr):
    c1_s4_sr.apply_photometric_component_spectral_response("dichroic")
    vec = c1_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_dichroic_c2(c2_s4_sr):
    c2_s4_sr.apply_photometric_component_spectral_response("dichroic")
    vec = c2_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_dichroic_c3(c3_s4_sr):
    c3_s4_sr.apply_photometric_component_spectral_response("dichroic")
    vec = c3_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_dichroic_c4(c4_s4_sr):
    c4_s4_sr.apply_photometric_component_spectral_response("dichroic")
    vec = c4_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_camera_abs(abs_s4_sr):
    abs_s4_sr.apply_photometric_component_spectral_response("camera")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_camera_c1(c1_s4_sr):
    c1_s4_sr.apply_photometric_component_spectral_response("camera")
    vec = c1_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_camera_c2(c2_s4_sr):
    c2_s4_sr.apply_photometric_component_spectral_response("camera")
    vec = c2_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_camera_c3(c3_s4_sr):
    c3_s4_sr.apply_photometric_component_spectral_response("camera")
    vec = c3_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_camera_c4(c4_s4_sr):
    c4_s4_sr.apply_photometric_component_spectral_response("camera")
    vec = c4_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_ccd_abs(abs_s4_sr):
    abs_s4_sr.apply_photometric_component_spectral_response("ccd")
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


def test_ccd_c1(c1_s4_sr):
    transmitance = [
        38.43,
        85.39,
        86.86,
        85.12,
        88.96,
        89.12,
        87.01,
        84.12,
        73.6,
        61.44,
        44.32,
        30.32,
        16.49,
        5.67,
        1.56,
        0.35,
    ]
    c1_s4_sr.apply_photometric_component_spectral_response("ccd")
    vec = c1_s4_sr.get_specific_flux()
    assert vec.all() == np.asarray(transmitance).all()


def test_ccd_c2(c2_s4_sr):
    transmitance = [
        38.98,
        50.04,
        75.82,
        95.8,
        95.77,
        92.86,
        88.07,
        81.16,
        68.9,
        56.67,
        42.7,
        27.75,
        14.9,
        6.24,
        1.25,
        0.27,
    ]
    c2_s4_sr.apply_photometric_component_spectral_response("ccd")
    vec = c2_s4_sr.get_specific_flux()
    assert vec.all() == np.asarray(transmitance).all()


def test_ccd_c3(c3_s4_sr):
    transmitance = [
        41.01,
        57.93,
        71.43,
        86.17,
        86.34,
        87.13,
        85.78,
        85.07,
        74.61,
        62.76,
        47.36,
        30.75,
        15.8,
        6.58,
        1.33,
        0.35,
    ]
    c3_s4_sr.apply_photometric_component_spectral_response("ccd")
    vec = c3_s4_sr.get_specific_flux()
    assert vec.all() == np.asarray(transmitance).all()


def test_ccd_c4(c4_s4_sr):
    transmitance = [
        71.33,
        79.5,
        88.78,
        96.43,
        89.44,
        88.56,
        82.91,
        69.85,
        63.45,
        57.65,
        40.06,
        26.09,
        13.91,
        6.65,
        1.4,
        0.34,
    ]
    c4_s4_sr.apply_photometric_component_spectral_response("ccd")
    vec = c4_s4_sr.get_specific_flux()
    assert vec.all() == np.asarray(transmitance).all()


# --------------------write specific_flux--------------------


def test_write_specific_flux(abs_s4_sr):
    specific_flux = np.asanyarray(
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
    )
    assert abs_s4_sr.get_specific_flux().all() == specific_flux.all()


# ---------------------- get_specific_flux -----------------------------


def test_get_specific_flux(abs_s4_sr):
    vec = abs_s4_sr.get_specific_flux()
    assert vec.all() == specific_flux.all()


# ----------------------- read_spreadsheet---------------------------


def test_read_spreadsheet_calibration_wheel(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/calibration_wheel.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_retarder(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/retarder.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_analyser(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/analyser.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_collimator(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/collimator.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 0/dichroic.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 0/camera.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 0/ccd.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_1(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 1/dichroic.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_1(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 1/camera.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_1(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 1/ccd.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 2/dichroic.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_2(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 2/camera.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_2(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 2/ccd.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_3(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 3/dichroic.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_3(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 3/camera.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_3(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 3/ccd.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_dichroic_4(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 4/dichroic.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_camera_4(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 4/camera.xlsx"
    abs_s4_sr._read_spreadsheet(file)


def test_read_spreadsheet_ccd_4(abs_s4_sr):
    file = "./SPARC4_Spectral_Response/Channel 4/ccd.xlsx"
    abs_s4_sr._read_spreadsheet(file)


# ----------------------- miscelaneous ----------------------------


def test_multiply_matrices(abs_s4_sr):
    a = np.ones((4, 4))
    abs_s4_sr._multiply_matrices(a, a)

    assert abs_s4_sr.specific_flux.all() == a.all()


def test_calculate_spline():
    transmitance = np.ones((1, n))[0]
    chc = Abstract_SPARC4_Spectral_Response()
    chc.write_specific_flux(specific_flux, wavelength_interval)
    new_transmitance = chc._calculate_spline(transmitance, wavelength_interval)
    assert new_transmitance.all() == transmitance.all()
