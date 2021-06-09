# -*- coding: utf-8 -*-
"""AIS strucute tests.

This script was created to test the structure of the classes as presented in
the UML diagram of the AIS project

Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Background_Image import Background_Image
from AIS.Channel_Creator import (
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)
from AIS.Header import Header
from AIS.Point_Spread_Function import Point_Spread_Function

import pytest

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 200,
}


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        star_magnitude=15,
        sky_magnitude=20,
        gaussian_std=3,
        ccd_operation_mode=dic,
        channel=1,
    )


@pytest.fixture
def chc1():
    return Concrete_Channel_1(ccd_temp=-70, sparc4_acquisition_mode="phot")


@pytest.fixture
def psf(chc1):
    return Point_Spread_Function(dic, 3, 3)


@pytest.fixture
def bgi(chc1):
    return Background_Image(dic, 3, 500)


# -------------------- Testing the AIS structure -----------------------------


def test_Channel_Creator(ais):
    var = 0
    if ais.CHC:
        var = 1
    assert var == 1


def test_Point_Spread_Function(ais):
    var = 0
    if ais.PSF:
        var = 1
    assert var == 1


def test_Background_Image(ais):
    var = 0
    if ais.BGI:
        var = 1
    assert var == 1


def test_HDR(ais):
    var = 0
    if ais.HDR:
        var = 1
    assert var == 1


def test_Spectrum_Calculation(ais):
    var = 0
    if ais.SC:
        var = 1
    assert var == 1


def test_Telescope_Spectral_Response(ais):
    var = 0
    if ais.TSR:
        var = 1
    assert var == 1


def test_Atmosphere_Spectral_Response(ais):
    var = 0
    if ais.ASR:
        var = 1
    assert var == 1
