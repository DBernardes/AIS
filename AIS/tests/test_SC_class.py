# -*- coding: utf-8 -*-
"""Spectrum calculation tests.

This script tests the operation of the spectrum calculation class.
"""

import numpy as np
import pytest
from AIS.Spectrum_Calculation import Spectrum_Calculation


@pytest.fixture
def sc():
    return Spectrum_Calculation(temperature=5700, l_init=350, l_final=1100, l_step=50)


def test_calculate_sky_specific_flux(sc):
    temp = [
        16681817027361.914,
        21112233469518.605,
        23675166849274.37,
        24574299823901.527,
        24243980202747.117,
        23117912376673.82,
        21542128053244.836,
        19762203592798.703,
        17938927532399.613,
        16170216117279.113,
        14510691893304.955,
        12986672615255.771,
        11606781118609.305,
        10369102934862.508,
        9265851307667.068,
    ]
    n = len(temp)
    specific_flux = np.zeros((4, n))
    specific_flux[0, :] = temp

    sky_specific_flux = sc.calculate_sky_specific_flux()
    for i in range(n):
        assert sky_specific_flux[1, i] == (specific_flux[1, i] * 0.1)


def test_calculate_star_specific_flux(sc):
    temp = [
        16681817027361.914,
        21112233469518.605,
        23675166849274.37,
        24574299823901.527,
        24243980202747.117,
        23117912376673.82,
        21542128053244.836,
        19762203592798.703,
        17938927532399.613,
        16170216117279.113,
        14510691893304.955,
        12986672615255.771,
        11606781118609.305,
        10369102934862.508,
        9265851307667.068,
    ]
    n = len(temp)
    specific_flux = np.zeros((4, n))
    specific_flux[0, :] = temp
    star_specific_flux = sc.calculate_star_specific_flux()
    for i in range(n):
        assert star_specific_flux[0, i] == specific_flux[0, i]
