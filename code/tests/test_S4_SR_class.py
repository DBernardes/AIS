# -*- coding: utf-8 -*-
"""SPARC4 spectrum response tests.

This script tests the operation of the SPARC4 spectrum response classes.
"""

from S4_SR import (Abstract_SPARC4_Spectral_Response,
                   Concrete_SPARC4_Spectral_Response_1,
                   Concrete_SPARC4_Spectral_Response_2,
                   Concrete_SPARC4_Spectral_Response_3,
                   Concrete_SPARC4_Spectral_Response_4)

import pytest


@pytest.fixture
def abs_s4_sr():
    return Abstract_SPARC4_Spectral_Response([1, 2, 3, 4, 5])


@pytest.fixture
def c1_s4_sr():
    return Concrete_SPARC4_Spectral_Response_1([1, 2, 3, 4, 5])


@pytest.fixture
def c2_s4_sr():
    return Concrete_SPARC4_Spectral_Response_2([1, 2, 3, 4, 5])


@pytest.fixture
def c3_s4_sr():
    return Concrete_SPARC4_Spectral_Response_3([1, 2, 3, 4, 5])


@pytest.fixture
def c4_s4_sr():
    return Concrete_SPARC4_Spectral_Response_4([1, 2, 3, 4, 5])


# -------------------- Initialize the class -----------------------

def test_spectrum_abs(abs_s4_sr):
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]


def test_spectrum_c1(c1_s4_sr):
    assert c1_s4_sr.spectrum == [1, 2, 3, 4, 5]


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
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- retarder  -----------------------


def test_retarder(abs_s4_sr):
    abs_s4_sr.retarder()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- analyzer -----------------------


def test_analyzer(abs_s4_sr):
    abs_s4_sr.analyzer()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- collimator -----------------------


def test_collimator(abs_s4_sr):
    abs_s4_sr.collimator()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- dichroic -----------------------


def test_dichroic(abs_s4_sr):
    abs_s4_sr.dichroic()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- camera -----------------------


def test_camera(abs_s4_sr):
    abs_s4_sr.camera()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- CCD -----------------------


def test_ccd(abs_s4_sr):
    abs_s4_sr.ccd()
    assert abs_s4_sr.spectrum == [1, 2, 3, 4, 5]

# -------------------- Integrate spectrum -----------------------


def test_integrate_spectrum(abs_s4_sr):
    flux = abs_s4_sr.integrate_spectrum()
    assert flux == sum([1, 2, 3, 4, 5])
