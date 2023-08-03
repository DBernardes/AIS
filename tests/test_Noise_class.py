# -*- coding: utf-8 -*-

"""
Tests of the Noise Class.

Created on Fri Apr 23 11:59:02 2021

@author: denis
"""

import pytest
import os
from AIS.Noise import Noise


@pytest.fixture
def noise():
    return Noise(1)


# ------------------------ Initialize the class --------------------------


def test_channel(noise):
    assert noise.channel == 1


def test_spreadsheet_path(noise):
    assert noise.spreadsheet_path == os.path.join(
        "AIS", "Noise", "spreadsheet", "Channel 1"
    )


# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, readout, preamp, binn, read_noise",
    [
        ("Conv", 1, 0.1, 1, 1, 8.78),
        ("Conv", 1, 0.1, 1, 2, 8.84),
        ("Conv", 1, 0.1, 2, 1, 3.46),
        ("Conv", 1, 0.1, 2, 2, 3.27),
        ("Conv", 1, 1, 1, 1, 6.67),
        ("Conv", 1, 1, 1, 2, 6.94),
        ("Conv", 1, 1, 2, 1, 4.76),
        ("Conv", 1, 1, 2, 2, 4.79),
    ],
)
def test_calc_read_noise_conv(
    noise, em_mode, em_gain, readout, preamp, binn, read_noise
):
    noise.em_mode = em_mode
    noise.em_gain = em_gain
    noise.readout = readout
    noise.preamp = preamp
    noise.binn = binn
    noise._calculate_read_noise_conventional_mode()
    assert round(noise.read_noise, 2) == read_noise


@pytest.mark.parametrize(
    "em_mode, em_gain, readout, preamp, binn, read_noise",
    [
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
def test_calc_read_noise_EM(noise, em_mode, em_gain, readout, preamp, binn, read_noise):
    noise.em_mode = em_mode
    noise.em_gain = em_gain
    noise.readout = readout
    noise.preamp = preamp
    noise.binn = binn
    noise._calculate_read_noise_em_mode()
    assert round(noise.read_noise, 2) == read_noise


@pytest.mark.parametrize(
    "em_mode, em_gain, readout, preamp, binn, read_noise",
    [
        ("EM", 2, 1, 1, 1, 24.64),
        ("Conv", 1, 0.1, 1, 1, 8.78),
    ],
)
def test_calc_read_noise(noise, em_mode, em_gain, readout, preamp, binn, read_noise):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "binn": binn,
        "preamp": preamp,
        "readout": readout,
    }
    rn = noise.calculate_read_noise(ccd_operation_mode)
    assert round(rn, 2) == read_noise


def test_calculate_dark_current_1(noise):
    dark_noise = noise.calculate_dark_current(-70)
    assert dark_noise == 5.8597559895090484e-05


def test_calculate_dark_current_2():
    noise = Noise(2)
    dark_noise = noise.calculate_dark_current(-70)
    assert dark_noise == 0.0001466809375420809


def test_calculate_dark_current_3():
    noise = Noise(3)
    dark_noise = noise.calculate_dark_current(-70)
    assert dark_noise == 8.688095666213275e-05


def test_calculate_dark_current_4():
    noise = Noise(4)
    dark_noise = noise.calculate_dark_current(-70)
    assert dark_noise == 0.00023133040352019222
