# -*- coding: utf-8 -*-

"""
SPARC4 spectral response tests.

This script tests the operation of the SPARC4 spectral response classes.
"""

import os

import numpy as np
import pytest
from AIS.SPARC4_Spectral_Response import (
    Abstract_SPARC4_Spectral_Response,
    Concrete_SPARC4_Spectral_Response_1,
    Concrete_SPARC4_Spectral_Response_2,
    Concrete_SPARC4_Spectral_Response_3,
    Concrete_SPARC4_Spectral_Response_4,
)

from .AIS_spectral_response_curves import (
    analyser_extra_ordinary_ray,
    analyser_ordinary_ray,
    calibration_wheel,
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
    dichroic_c0,
    dichroic_c1,
    dichroic_c2,
    dichroic_c3,
    dichroic_c4,
    retarder,
    star_specific_flux,
    wavelength_interval,
    wavelength_interval_len,
)

# -------------------------------------------------------------------------------------------------------------


@pytest.fixture
def abs_s4_sr():
    chc = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    return chc


@pytest.fixture
def c1_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_1(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    return chc


@pytest.fixture
def c2_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_2(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    return chc


@pytest.fixture
def c3_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_3(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    return chc


@pytest.fixture
def c4_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_4(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    return chc


def multiply_matrices(matrix, star_specific_flux):
    for i in range(len(star_specific_flux[0])):
        star_specific_flux[:, i] = np.dot(matrix, star_specific_flux[:, i])
    return star_specific_flux


# -------------------- Initialize the class -----------------------


def test_specific_ordinary_ray_abs(abs_s4_sr):
    abs_specific_ordinary_ray = abs_s4_sr.get_specific_ordinary_ray()
    assert np.allclose(abs_specific_ordinary_ray, star_specific_flux)


def test_specific_ordinary_ray_c1(c1_s4_sr):
    c1_specific_ordinary_ray = c1_s4_sr.get_specific_ordinary_ray()
    assert np.allclose(c1_specific_ordinary_ray, star_specific_flux)


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


ordinary_ray = multiply_matrices(analyser_ordinary_ray, star_specific_flux.copy())
extra_ordinary_ray = multiply_matrices(
    analyser_extra_ordinary_ray, star_specific_flux.copy()
)


def test_calibration_wheel(abs_s4_sr):
    abs_s4_sr.apply_calibration_wheel()
    new_star_specific_flux = multiply_matrices(calibration_wheel, star_specific_flux)
    assert np.allclose(abs_s4_sr.get_specific_ordinary_ray(), new_star_specific_flux)


def test_retarder(abs_s4_sr):
    abs_s4_sr.apply_retarder()
    new_star_specific_flux = multiply_matrices(retarder, star_specific_flux)
    assert np.allclose(abs_s4_sr.get_specific_ordinary_ray(), new_star_specific_flux)


def test_analyzer(abs_s4_sr):
    abs_s4_sr.apply_analyser()
    assert np.allclose(abs_s4_sr.specific_ordinary_ray, ordinary_ray)
    assert np.allclose(abs_s4_sr.specific_extra_ordinary_ray, extra_ordinary_ray)


def test_collimator(abs_s4_sr):
    abs_s4_sr.apply_analyser()
    abs_s4_sr.apply_collimator()
    new_ordinary_ray = np.multiply(ordinary_ray.copy(), colimator_transmitance)
    new_extra_ordinary_ray = np.multiply(
        extra_ordinary_ray.copy(), colimator_transmitance
    )
    assert np.allclose(abs_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(abs_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_dichroic_abs(abs_s4_sr):
    abs_s4_sr.apply_analyser()
    abs_s4_sr.apply_dichroic()
    new_ordinary_ray = np.multiply(ordinary_ray.copy(), dichroic_c0)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), dichroic_c0)

    assert np.allclose(abs_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(abs_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_dichroic_c1(c1_s4_sr):
    c1_s4_sr.apply_analyser()
    c1_s4_sr.apply_dichroic()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), dichroic_c1)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), dichroic_c1)

    assert np.allclose(c1_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c1_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_dichroic_c2(c2_s4_sr):
    c2_s4_sr.apply_analyser()
    c2_s4_sr.apply_dichroic()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), dichroic_c2)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), dichroic_c2)

    assert np.allclose(c2_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c2_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_dichroic_c3(c3_s4_sr):
    c3_s4_sr.apply_analyser()
    c3_s4_sr.apply_dichroic()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), dichroic_c3)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), dichroic_c3)

    assert np.allclose(c3_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c3_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_dichroic_c4(c4_s4_sr):
    c4_s4_sr.apply_analyser()
    c4_s4_sr.apply_dichroic()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), dichroic_c4)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), dichroic_c4)

    assert np.allclose(c4_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c4_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_camera_abs(abs_s4_sr):
    abs_s4_sr.apply_analyser()
    abs_s4_sr.apply_camera()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), camera_c0)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), camera_c0)

    assert np.allclose(abs_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(abs_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_camera_c1(c1_s4_sr):
    c1_s4_sr.apply_analyser()
    c1_s4_sr.apply_camera()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), camera_c1)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), camera_c1)

    assert np.allclose(c1_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c1_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_camera_c2(c2_s4_sr):
    c2_s4_sr.apply_analyser()
    c2_s4_sr.apply_camera()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), camera_c2)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), camera_c2)

    assert np.allclose(c2_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c2_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_camera_c3(c3_s4_sr):
    c3_s4_sr.apply_analyser()
    c3_s4_sr.apply_camera()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), camera_c3)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), camera_c3)

    assert np.allclose(c3_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c3_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_camera_c4(c4_s4_sr):
    c4_s4_sr.apply_analyser()
    c4_s4_sr.apply_camera()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), camera_c4)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), camera_c4)

    assert np.allclose(c4_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c4_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_ccd_abs(abs_s4_sr):
    abs_s4_sr.apply_analyser()
    abs_s4_sr.apply_ccd()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), ccd_transmitance_c0)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), ccd_transmitance_c0)

    assert np.allclose(abs_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(abs_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_ccd_c1(c1_s4_sr):
    c1_s4_sr.apply_analyser()
    c1_s4_sr.apply_ccd()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), ccd_transmitance_c1)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), ccd_transmitance_c1)

    assert np.allclose(c1_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c1_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_ccd_c2(c2_s4_sr):
    c2_s4_sr.apply_analyser()
    c2_s4_sr.apply_ccd()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), ccd_transmitance_c2)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), ccd_transmitance_c2)

    assert np.allclose(c2_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c2_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_ccd_c3(c3_s4_sr):
    c3_s4_sr.apply_analyser()
    c3_s4_sr.apply_ccd()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), ccd_transmitance_c3)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), ccd_transmitance_c3)

    assert np.allclose(c3_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c3_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


def test_ccd_c4(c4_s4_sr):
    c4_s4_sr.apply_analyser()
    c4_s4_sr.apply_ccd()

    new_ordinary_ray = np.multiply(ordinary_ray.copy(), ccd_transmitance_c4)
    new_extra_ordinary_ray = np.multiply(extra_ordinary_ray.copy(), ccd_transmitance_c4)

    assert np.allclose(c4_s4_sr.specific_ordinary_ray, new_ordinary_ray)
    assert np.allclose(c4_s4_sr.specific_extra_ordinary_ray, new_extra_ordinary_ray)


# --------------------write star_specific_flux--------------------


def test_write_specific_flux():
    star_specific_flux = np.asanyarray(
        [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
    )
    wavelength_interval = range(400, 1150, 50)
    s4_sr = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    s4_sr.write_specific_flux(star_specific_flux)
    assert np.allclose(s4_sr.specific_ordinary_ray, star_specific_flux)


# ---------------------- get specific ordinary and extraordinary rays -----------------------------


def test_get_specific_ordinary_ray(abs_s4_sr):
    abs_ordinary_ray = abs_s4_sr.get_specific_ordinary_ray()
    assert np.allclose(abs_ordinary_ray, star_specific_flux)


def test_get_specific_extra_ordinary_ray(abs_s4_sr):
    assert np.allclose(
        abs_s4_sr.get_specific_extra_ordinary_ray(),
        np.zeros((4, wavelength_interval_len)),
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
    star_specific_flux = abs_s4_sr._multiply_matrices(a, a)
    boolean_test = star_specific_flux == a
    assert boolean_test.all()


def test_calculate_spline():
    transmitance = np.ones((1, wavelength_interval_len))[0]
    chc = Abstract_SPARC4_Spectral_Response(wavelength_interval)
    chc.write_specific_flux(star_specific_flux)
    new_transmitance = chc._calculate_spline(transmitance, wavelength_interval)
    assert np.allclose(new_transmitance, transmitance)
