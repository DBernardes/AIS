# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 13:37:22 2021

@author: denis
"""

# -------------------- Testing the AIS class --------------------------------
from AIS import Artificial_Images_Simulator
import pytest


@pytest.fixture
def ais():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    return Artificial_Images_Simulator(100.0, 10.0, 3, dic)


dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1, 'hss': 1, 'bin': 1, 't_exp': 1}


# ---------------------_calc_dark_current--------------------------------


def test_dark_current_9914(ais):
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(5.85975599e-5, 7)


def test_dark_current_9915(ais):
    ais.serial_number = 9915
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(1.466809375e-4, 7)


def test_dark_current_9916(ais):
    ais.serial_number = 9916
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(8.688095666e-5, 7)


def test_dark_current_9917(ais):
    ais.serial_number = 9917
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(2.313304035e-4, 7)

# -----------------------------_calc_read_noise--------------------------


@pytest.mark.parametrize(
    'em_mode, em_gain, hss, preamp, binn, read_noise',
    [(0, 2, 0.1, 1, 1, 8.78),
     (0, 2, 0.1, 1, 2, 8.84),
     (0, 2, 0.1, 2, 1, 3.46),
     (0, 2, 0.1, 2, 2, 3.27),
     (0, 2, 1, 1, 1, 6.67),
     (0, 2, 1, 1, 2, 6.94),
     (0, 2, 1, 2, 1, 4.76),
     (0, 2, 1, 2, 2, 4.79),
     (1, 2, 1, 1, 1, 24.64),
     (1, 2, 1, 1, 2, 33.76),
     (1, 2, 1, 2, 1, 12.05),
     (1, 2, 1, 2, 2, 14.55),
     (1, 2, 10, 1, 1, 83.68),
     (1, 2, 10, 1, 2, 82.93),
     (1, 2, 10, 2, 1, 41.71),
     (1, 2, 10, 2, 2, 41.82),
     (1, 2, 20, 1, 1, 160.06),
     (1, 2, 20, 1, 2, 161.98),
     (1, 2, 20, 2, 1, 66.01),
     (1, 2, 20, 2, 2, 72.71),
     (1, 2, 30, 1, 1, 262.01),
     (1, 2, 30, 1, 2, 273.19),
     (1, 2, 30, 2, 1, 169.25),
     (1, 2, 30, 2, 2, 143.59),
     ]
)
def test_calc_read_noise(ais, em_mode, em_gain, hss, preamp, binn, read_noise):
    ais.em_mode = em_mode
    ais.em_gain = em_gain
    ais.hss = hss
    ais.preamp = preamp
    ais.preamp = preamp
    ais.bin = binn
    ais._calc_read_noise()
    assert round(ais.read_noise, 2) == read_noise
