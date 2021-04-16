# -*- coding: utf-8 -*-
"""AIS strucute tests.

This script was created to test the structure of the classes as presented in
the UML diagram of the AIS project

Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


from AIS import Artificial_Images_Simulator
from PSF import Point_Spread_Function
from BGI import Background_Image
from HDR import Header
from CHC import (Concrete_Channel_1,
                 Concrete_Channel_2,
                 Concrete_Channel_3,
                 Concrete_Channel_4)
import pytest


dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
       'hss': 1, 'bin': 1, 't_exp': 1, 'ccd_temp': -70}


@pytest.fixture
def ais():
    return Artificial_Images_Simulator(star_flux=100.0,
                                       sky_flux=10.0,
                                       gaussian_stddev=3,
                                       ccd_operation_mode=dic)


@pytest.fixture
def chc1():
    return Concrete_Channel_1()


@pytest.fixture
def psf(chc1):
    return Point_Spread_Function(chc1)


@pytest.fixture
def bgi(chc1):
    return Background_Image(chc1)


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


# -------------------- Testing the PSF structure -----------------------------


def test_Point_Spread_Function_Channel_Creator(psf):
    var = 0
    if psf.CHC:
        var = 1
    assert var == 1


def test_Point_Spread_Function_Flux_Calculation(psf):
    var = 0
    if psf.FC:
        var = 1
    assert var == 1


def test_Point_Spread_Function_Telescope_Spectral_Response(psf):
    var = 0
    if psf.TSR:
        var = 1
    assert var == 1


def test_Point_Spread_Function_Atmosphere_Spectral_Response(psf):
    var = 0
    if psf.ASR:
        var = 1
    assert var == 1

# -------------------- Testing the BGI structure -----------------------------


def test_Background_Image_Channel_Creator(bgi):
    var = 0
    if bgi.CHC:
        var = 1
    assert var == 1


def test_Background_Image_Flux_Calculation(bgi):
    var = 0
    if bgi.FC:
        var = 1
    assert var == 1


def test_Background_Imagen_Telescope_Spectral_Response(bgi):
    var = 0
    if bgi.TSR:
        var = 1
    assert var == 1


def test_Background_Image_Atmosphere_Spectral_Response(bgi):
    var = 0
    if bgi.ASR:
        var = 1
    assert var == 1
