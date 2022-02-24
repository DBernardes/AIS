# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""

import numpy as np
import pytest
from AIS.Channel_Creator import (
    Abstract_Channel_Creator,
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)
from AIS.Spectrum_Calculation import Spectrum_Calculation

from .AIS_spectral_response_curves import (
    POLARIZER_90_MATRIX,
    POLARIZER_MATRIX,
    THETA_POL,
    analyser_transmitance,
    calculate_retarder_matrix,
    camera_c1,
    camera_c2,
    camera_c3,
    camera_c4,
    ccd_operation_mode,
    ccd_temp,
    ccd_transmitance_c1,
    ccd_transmitance_c2,
    ccd_transmitance_c3,
    ccd_transmitance_c4,
    collimator_transmitance,
    depolarizer_transmitance,
    dichroic_c1,
    dichroic_c2,
    dichroic_c3,
    dichroic_c4,
    polarizer_transmitance,
    retardance_half,
    retardance_quarter,
    retarder_transmitance,
    sparc4_operation_mode,
    star_specific_photons_per_second,
    wavelength_interval,
)
from .test_AIS_operation import multiply_matrices

# --------------------------------------------------------------------------------------------------------------------


@pytest.fixture
def abs_chc():
    return Abstract_Channel_Creator(
        wavelength_interval=wavelength_interval,
        ccd_operation_mode=ccd_operation_mode,
    )


@pytest.fixture
def chc1():
    return Concrete_Channel_1(
        wavelength_interval=wavelength_interval,
        ccd_operation_mode=ccd_operation_mode,
    )


@pytest.fixture
def chc2():
    return Concrete_Channel_2(
        wavelength_interval=wavelength_interval,
        ccd_operation_mode=ccd_operation_mode,
    )


@pytest.fixture
def chc3():
    return Concrete_Channel_3(
        wavelength_interval=wavelength_interval,
        ccd_operation_mode=ccd_operation_mode,
    )


@pytest.fixture
def chc4():
    return Concrete_Channel_4(
        wavelength_interval=wavelength_interval,
        ccd_operation_mode=ccd_operation_mode,
    )


# ------------------------ Initialize the class --------------------------


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


# ----------------------- Calculate read noise  -------------------------


def test_calculate_read_noise_1(chc1):
    chc1.calculate_read_noise()
    assert chc1.read_noise == 6.67


def test_calculate_read_noise_2(chc2):
    chc2.calculate_read_noise()
    assert chc2.read_noise == 6.67


def test_calculate_read_noise_3(chc3):
    chc3.calculate_read_noise()
    assert chc3.read_noise == 6.67


def test_calculate_read_noise_4(chc4):
    chc4.calculate_read_noise()
    assert chc4.read_noise == 6.67


# ----------------------- Test CCD gain -------------------------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, ccd_gain",
    [
        ("Conv", 1, 0.1, 1, 1, 3.35),
        ("Conv", 1, 0.1, 2, 1, 0.80),
        ("Conv", 1, 1, 1, 1, 3.37),
        ("Conv", 1, 1, 2, 1, 0.80),
        ("EM", 2, 1, 1, 1, 15.90),
        ("EM", 2, 1, 2, 1, 3.88),
        ("EM", 2, 10, 1, 1, 16.00),
        ("EM", 2, 10, 2, 1, 3.96),
        ("EM", 2, 20, 1, 1, 16.40),
        ("EM", 2, 20, 2, 1, 4.39),
        ("EM", 2, 30, 1, 1, 17.20),
        ("EM", 2, 30, 2, 1, 5.27),
    ],
)
def test_configure_gain_ch1(chc1, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    chc1.ccd_operation_mode = ccd_operation_mode
    ccd_gain = chc1.get_ccd_gain()
    assert ccd_gain == ccd_gain


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, ccd_gain",
    [
        ("Conv", 1, 0.1, 1, 1, 3.35),
        ("Conv", 1, 0.1, 2, 1, 0.80),
        ("Conv", 1, 1, 1, 1, 3.37),
        ("Conv", 1, 1, 2, 1, 0.80),
        ("EM", 2, 1, 1, 1, 15.90),
        ("EM", 2, 1, 2, 1, 3.88),
        ("EM", 2, 10, 1, 1, 16.00),
        ("EM", 2, 10, 2, 1, 3.96),
        ("EM", 2, 20, 1, 1, 16.40),
        ("EM", 2, 20, 2, 1, 4.39),
        ("EM", 2, 30, 1, 1, 17.20),
        ("EM", 2, 30, 2, 1, 5.27),
    ],
)
def test_configure_gain_ch2(chc2, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    chc2.ccd_operation_mode = ccd_operation_mode
    ccd_gain = chc2.get_ccd_gain()
    assert ccd_gain == ccd_gain


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, ccd_gain",
    [
        ("Conv", 1, 0.1, 1, 1, 3.35),
        ("Conv", 1, 0.1, 2, 1, 0.80),
        ("Conv", 1, 1, 1, 1, 3.37),
        ("Conv", 1, 1, 2, 1, 0.80),
        ("EM", 2, 1, 1, 1, 15.90),
        ("EM", 2, 1, 2, 1, 3.88),
        ("EM", 2, 10, 1, 1, 16.00),
        ("EM", 2, 10, 2, 1, 3.96),
        ("EM", 2, 20, 1, 1, 16.40),
        ("EM", 2, 20, 2, 1, 4.39),
        ("EM", 2, 30, 1, 1, 17.20),
        ("EM", 2, 30, 2, 1, 5.27),
    ],
)
def test_configure_gain_ch3(chc3, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    chc3.ccd_operation_mode = ccd_operation_mode
    ccd_gain = chc3.get_ccd_gain()
    assert ccd_gain == ccd_gain


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, ccd_gain",
    [
        ("Conv", 1, 0.1, 1, 1, 3.35),
        ("Conv", 1, 0.1, 2, 1, 0.80),
        ("Conv", 1, 1, 1, 1, 3.37),
        ("Conv", 1, 1, 2, 1, 0.80),
        ("EM", 2, 1, 1, 1, 15.90),
        ("EM", 2, 1, 2, 1, 3.88),
        ("EM", 2, 10, 1, 1, 16.00),
        ("EM", 2, 10, 2, 1, 3.96),
        ("EM", 2, 20, 1, 1, 16.40),
        ("EM", 2, 20, 2, 1, 4.39),
        ("EM", 2, 30, 1, 1, 17.20),
        ("EM", 2, 30, 2, 1, 5.27),
    ],
)
def test_configure_gain_ch4(chc4, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    chc4.ccd_operation_mode = ccd_operation_mode
    ccd_gain = chc4.get_ccd_gain()
    assert ccd_gain == ccd_gain


# ------------------- Teste apply_sparc4_spectral_response - quarter waveplate  -----------------------------


def test_apply_sparc4_spectral_response_retarder_quarter():
    # Applying retarder
    star_specific_photons_per_second_c1 = []
    for array in star_specific_photons_per_second.copy():
        for idx, value in enumerate(array[0]):
            retarder_matrix = calculate_retarder_matrix(
                retardance_quarter[idx], THETA_POL
            )
            array[:, idx] = np.dot(retarder_matrix, array[:, idx])
        star_specific_photons_per_second_c1.append(array)
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        retarder_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying analyser
    ordinary_ray = multiply_matrices(
        POLARIZER_MATRIX, star_specific_photons_per_second_c1.copy()
    )
    ordinary_ray[0][0] = np.multiply(analyser_transmitance, ordinary_ray[0][0])
    extra_ordinary_ray = multiply_matrices(
        POLARIZER_90_MATRIX, star_specific_photons_per_second_c1
    )
    extra_ordinary_ray[0][0] = np.multiply(
        analyser_transmitance, extra_ordinary_ray[0][0]
    )
    star_specific_photons_per_second_c1 = [ordinary_ray[0], extra_ordinary_ray[0]]
    # Applying collimator
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        collimator_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying camera
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        camera_c1, star_specific_photons_per_second_c1[0][0]
    )
    # Applying ccd
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
    )

    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
    }
    chc1_pol = Concrete_Channel_1(wavelength_interval, ccd_operation_mode)
    new_star_specific_photons_per_second = chc1_pol.apply_sparc4_spectral_response(
        star_specific_photons_per_second.copy(), sparc4_operation_mode
    )
    assert np.allclose(
        new_star_specific_photons_per_second, star_specific_photons_per_second_c1
    )


# ------------------- Teste apply_sparc4_spectral_response - half waveplate  -----------------------------


def test_apply_sparc4_spectral_response_retarder_half():
    # Applying retarder
    star_specific_photons_per_second_c1 = []
    for array in star_specific_photons_per_second.copy():
        for idx, value in enumerate(array[0]):
            retarder_matrix = calculate_retarder_matrix(retardance_half[idx], THETA_POL)
            array[:, idx] = np.dot(retarder_matrix, array[:, idx])
        star_specific_photons_per_second_c1.append(array)
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        retarder_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying analyser
    ordinary_ray = multiply_matrices(
        POLARIZER_MATRIX, star_specific_photons_per_second_c1.copy()
    )
    ordinary_ray[0][0] = np.multiply(analyser_transmitance, ordinary_ray[0][0])
    extra_ordinary_ray = multiply_matrices(
        POLARIZER_90_MATRIX, star_specific_photons_per_second_c1
    )
    extra_ordinary_ray[0][0] = np.multiply(
        analyser_transmitance, extra_ordinary_ray[0][0]
    )
    star_specific_photons_per_second_c1 = [ordinary_ray[0], extra_ordinary_ray[0]]
    # Applying collimator
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        collimator_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying camera
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        camera_c1, star_specific_photons_per_second_c1[0][0]
    )
    # Applying ccd
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
    )

    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "half",
    }
    chc1_pol = Concrete_Channel_1(wavelength_interval, ccd_operation_mode)
    new_star_specific_photons_per_second = chc1_pol.apply_sparc4_spectral_response(
        star_specific_photons_per_second.copy(), sparc4_operation_mode
    )
    assert np.allclose(
        new_star_specific_photons_per_second, star_specific_photons_per_second_c1
    )


# ------------------- Teste apply_sparc4_spectral_response - polarizer + quarter waveplate  -----------------------------


def test_apply_sparc4_spectral_response_retarder_quarter_polarizer():
    # Apply polarizer
    star_specific_photons_per_second_c1 = multiply_matrices(
        POLARIZER_MATRIX, star_specific_photons_per_second.copy()
    )
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        polarizer_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying retarder
    temp = []
    for array in star_specific_photons_per_second_c1:
        for idx, value in enumerate(array[0]):
            retarder_matrix = calculate_retarder_matrix(
                retardance_quarter[idx], THETA_POL
            )
            array[:, idx] = np.dot(retarder_matrix, array[:, idx])
        temp.append(array)
    star_specific_photons_per_second_c1 = temp
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        retarder_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying analyser
    ordinary_ray = multiply_matrices(
        POLARIZER_MATRIX, star_specific_photons_per_second_c1.copy()
    )
    ordinary_ray[0][0] = np.multiply(analyser_transmitance, ordinary_ray[0][0])
    extra_ordinary_ray = multiply_matrices(
        POLARIZER_90_MATRIX, star_specific_photons_per_second_c1
    )
    extra_ordinary_ray[0][0] = np.multiply(
        analyser_transmitance, extra_ordinary_ray[0][0]
    )
    star_specific_photons_per_second_c1 = [ordinary_ray[0], extra_ordinary_ray[0]]
    # Applying collimator
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        collimator_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying camera
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        camera_c1, star_specific_photons_per_second_c1[0][0]
    )
    # Applying ccd
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
    )

    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "polarizer",
        "retarder": "quarter",
    }
    chc1_pol = Concrete_Channel_1(wavelength_interval, ccd_operation_mode)
    new_star_specific_photons_per_second = chc1_pol.apply_sparc4_spectral_response(
        star_specific_photons_per_second.copy(), sparc4_operation_mode
    )
    assert np.allclose(
        new_star_specific_photons_per_second, star_specific_photons_per_second_c1
    )


# ------------------- Teste apply_sparc4_spectral_response - photometric  -----------------------------


def test_apply_sparc4_spectral_response_photometric():

    star_specific_photons_per_second_c1 = star_specific_photons_per_second.copy()
    # Applying collimator
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        collimator_transmitance, star_specific_photons_per_second_c1[0][0]
    )
    # Applying camera
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        camera_c1, star_specific_photons_per_second_c1[0][0]
    )
    # Applying ccd
    star_specific_photons_per_second_c1[0][0] = np.multiply(
        ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
    )

    sparc4_operation_mode = {
        "acquisition_mode": "photometric",
    }
    chc1_pol = Concrete_Channel_1(wavelength_interval, ccd_operation_mode)
    new_star_specific_photons_per_second = chc1_pol.apply_sparc4_spectral_response(
        star_specific_photons_per_second.copy(), sparc4_operation_mode
    )
    assert np.allclose(
        new_star_specific_photons_per_second, star_specific_photons_per_second_c1
    )
