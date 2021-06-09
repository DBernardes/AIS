# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""


from code.CHC.CHC import (
    Abstract_Channel_Creator,
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)

import pytest


@pytest.fixture
def abs_chc():
    return Abstract_Channel_Creator(ccd_temp=-70, sparc4_operation_mode="phot")


@pytest.fixture
def chc1():
    return Concrete_Channel_1(ccd_temp=-70, sparc4_operation_mode="phot")


@pytest.fixture
def chc2():
    return Concrete_Channel_2(ccd_temp=-70, sparc4_operation_mode="phot")


@pytest.fixture
def chc3():
    return Concrete_Channel_3(ccd_temp=-70, sparc4_operation_mode="phot")


@pytest.fixture
def chc4():
    return Concrete_Channel_4(ccd_temp=-70, sparc4_operation_mode="phot")


# ------------------------ Initialize the class --------------------------


def test_ccd_temp(abs_chc):
    assert abs_chc.ccd_temp == -70


def test_sparc_acquisition_mode(abs_chc):
    assert abs_chc.sparc4_operation_mode == "phot"


def test_channel_ID_abs(abs_chc):
    assert abs_chc._CHANNEL_ID == 0


def test_serial_number_abs(abs_chc):
    assert abs_chc._SERIAL_NUMBER == 0


def test_channel_ID_1(chc1):
    assert chc1._CHANNEL_ID == 1


def test_serial_number_1(chc1):
    assert chc1._SERIAL_NUMBER == 9914


def test_channel_ID_2(chc2):
    assert chc2._CHANNEL_ID == 2


def test_serial_number_2(chc2):
    assert chc2._SERIAL_NUMBER == 9915


def test_channel_ID_3(chc3):
    assert chc3._CHANNEL_ID == 3


def test_serial_number_3(chc3):
    assert chc3._SERIAL_NUMBER == 9916


def test_channel_ID_4(chc4):
    assert chc4._CHANNEL_ID == 4


def test_serial_number_4(chc4):
    assert chc4._SERIAL_NUMBER == 9917


# ----------------------- Calculate dark current  -------------------------


def test_calculate_dark_current_1(chc1):
    chc1.calculate_dark_current()
    assert round(chc1.dark_current, 7) == 5.86e-5


def test_calculate_dark_current_2(chc2):
    chc2.calculate_dark_current()
    assert chc2.dark_current, 7 == 0.0001467


def test_calculate_dark_current_3(chc3):
    chc3.calculate_dark_current()
    assert round(chc3.dark_current, 7) == 8.69e-05


def test_calculate_dark_current_4(chc4):
    chc4.calculate_dark_current()
    assert chc4.dark_current, 7 == 0.0002313


# # -------------------------Calculate Read Noise -------------------------

# def test_calculate_read_noise(bgi):
#     bgi._calculate_read_noise(dic)
#     assert bgi.read_noise == 6.67
