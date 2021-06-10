# -*- coding: utf-8 -*-
"""SPARC4 spectrum response tests.

This script tests the operation of the SPARC4 spectrum response classes.
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

spectrum = np.ones((1, 4))[0]
n = len(spectrum)


@pytest.fixture
def abs_s4_sr():
    chc = Abstract_SPARC4_Spectral_Response()
    chc.write_spectrum(spectrum)
    return chc


@pytest.fixture
def c1_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_1()
    chc.write_spectrum(spectrum)
    return chc


@pytest.fixture
def c2_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_2()
    chc.write_spectrum(spectrum)
    return chc


@pytest.fixture
def c3_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_3()
    chc.write_spectrum(spectrum)
    return chc


@pytest.fixture
def c4_s4_sr():
    chc = Concrete_SPARC4_Spectral_Response_4()
    chc.write_spectrum(spectrum)
    return chc


# -------------------- Initialize the class -----------------------


def test_spectrum_abs(abs_s4_sr):
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_spectrum_c1(c1_s4_sr):
    vec = c1_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


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


# -------------------- calibration wheel  -----------------------


def test_calibration_wheel(abs_s4_sr):
    abs_s4_sr.calibration_wheel()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- retarder  -----------------------


def test_retarder(abs_s4_sr):
    abs_s4_sr.retarder()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- analyzer -----------------------


def test_analyzer(abs_s4_sr):
    abs_s4_sr.analyzer()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- collimator -----------------------


def test_collimator(abs_s4_sr):
    abs_s4_sr.collimator()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- dichroic -----------------------


def test_dichroic_abs(abs_s4_sr):
    abs_s4_sr.dichroic()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_dichroic_c1(c1_s4_sr):
    c1_s4_sr.dichroic()
    vec = c1_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_dichroic_c2(c2_s4_sr):
    c2_s4_sr.dichroic()
    vec = c2_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_dichroic_c3(c3_s4_sr):
    c3_s4_sr.dichroic()
    vec = c3_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_dichroic_c4(c4_s4_sr):
    c4_s4_sr.dichroic()
    vec = c4_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- camera -----------------------


def test_camera_abs(abs_s4_sr):
    abs_s4_sr.camera()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_camera_c1(c1_s4_sr):
    c1_s4_sr.camera()
    vec = c1_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_camera_c2(c2_s4_sr):
    c2_s4_sr.camera()
    vec = c2_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_camera_c3(c3_s4_sr):
    c3_s4_sr.camera()
    vec = c3_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_camera_c4(c4_s4_sr):
    c4_s4_sr.camera()
    vec = c4_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# -------------------- CCD -----------------------


def test_ccd_abs(abs_s4_sr):
    abs_s4_sr.ccd()
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_ccd_c1(c1_s4_sr):
    c1_s4_sr.ccd()
    vec = c1_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_ccd_c2(c2_s4_sr):
    c2_s4_sr.ccd()
    vec = c2_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_ccd_c3(c3_s4_sr):
    c3_s4_sr.ccd()
    vec = c3_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


def test_ccd_c4(c4_s4_sr):
    c4_s4_sr.ccd()
    vec = c4_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


# --------------------write spectrum--------------------


def test_write_spectrum(abs_s4_sr):
    spectrum = [[1, 1, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    new_spectrum = [list(line) for line in abs_s4_sr.get_spectrum()]
    assert new_spectrum == spectrum


# ---------------------- get_spectrum -----------------------------


def test_get_spectrum(abs_s4_sr):
    vec = abs_s4_sr.get_spectrum()[0, :]
    for i in range(n):
        assert vec[i] == spectrum[i]


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
