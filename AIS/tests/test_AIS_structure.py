# -*- coding: utf-8 -*-

"""AIS strucute tests.

This script was created to test the structure of the classes as presented in
the UML diagram of the AIS project

Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

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
        ccd_operation_mode=dic,
        channel=1,
    )


# -------------------- Testing the AIS structure -----------------------------


def test_Channel_Creator(ais):
    var = 0
    if ais.chc:
        var = 1
    assert var == 1


def test_Point_Spread_Function(ais):
    var = 0
    if ais.psf:
        var = 1
    assert var == 1


def test_Background_Image(ais):
    var = 0
    if ais.bgi:
        var = 1
    assert var == 1


def test_HDR(ais):
    var = 0
    if ais.hdr:
        var = 1
    assert var == 1


def test_Spectrum_Calculation(ais):
    var = 0
    if ais.sc:
        var = 1
    assert var == 1


def test_Telescope_Spectral_Response(ais):
    var = 0
    if ais.tsr:
        var = 1
    assert var == 1


def test_Atmosphere_Spectral_Response(ais):
    var = 0
    if ais.asr:
        var = 1
    assert var == 1
