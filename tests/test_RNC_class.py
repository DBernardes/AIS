# -*- coding: utf-8 -*-

"""
Test of the Read Noise Calculation Class.

Created on Fri Apr 23 11:59:02 2021

@author: denis
"""

import pytest
from AIS.Read_Noise_Calculation import Read_Noise_Calculation

dic = {"em_mode": "Conv", "em_gain": 1, "binn": 1, "preamp": 1, "hss": 1}


@pytest.fixture
def rnc():
    return Read_Noise_Calculation(dic, 1)


# ------------------------ Initialize the class --------------------------


def test_em_mode(rnc):
    assert rnc.em_mode == "Conv"


def test_em_gain(rnc):
    assert rnc.em_gain == 1


def test_binn(rnc):
    assert rnc.binn == 1


def test_preamp(rnc):
    assert rnc.preamp == 1


def test_hss(rnc):
    assert rnc.hss == 1


def test_directory(rnc):
    assert rnc.directory == "Channel 1"


# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, read_noise",
    [
        ("Conv", 1, 0.1, 1, 1, 8.78),
        ("Conv", 1, 0.1, 1, 2, 8.84),
        ("Conv", 1, 0.1, 2, 1, 3.46),
        ("Conv", 1, 0.1, 2, 2, 3.27),
        ("Conv", 1, 1, 1, 1, 6.67),
        ("Conv", 1, 1, 1, 2, 6.94),
        ("Conv", 1, 1, 2, 1, 4.76),
        ("Conv", 1, 1, 2, 2, 4.79),
        ("EM", 2, 1, 1, 1, 24.64),
        ("EM", 2, 1, 1, 2, 33.76),
        ("EM", 2, 1, 2, 1, 12.05),
        ("EM", 2, 1, 2, 2, 14.55),
        ("EM", 2, 10, 1, 1, 83.68),
        ("EM", 2, 10, 1, 2, 82.93),
        ("EM", 2, 10, 2, 1, 41.71),
        ("EM", 2, 10, 2, 2, 41.82),
        ("EM", 2, 20, 1, 1, 160.06),
        ("EM", 2, 20, 1, 2, 161.98),
        ("EM", 2, 20, 2, 1, 66.01),
        ("EM", 2, 20, 2, 2, 72.71),
        ("EM", 2, 30, 1, 1, 262.01),
        ("EM", 2, 30, 1, 2, 273.19),
        ("EM", 2, 30, 2, 1, 169.25),
        ("EM", 2, 30, 2, 2, 143.59),
    ],
)
def test_calc_read_noise(rnc, em_mode, em_gain, hss, preamp, binn, read_noise):
    rnc.em_mode = em_mode
    rnc.em_gain = em_gain
    rnc.hss = hss
    rnc.preamp = preamp
    rnc.binn = binn
    rn = rnc.calculate_read_noise()
    assert round(rn, 2) == read_noise


def test_get_operation_mode(rnc):
    rnc.get_operation_mode()
