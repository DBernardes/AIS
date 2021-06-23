# -*- coding: utf-8 -*-
"""Flux Calculation class tests.

This script tests the operation of the Background Image Class.

Created on Thu Apr 22 13:44:35 2021

@author: denis
"""

import os

import numpy as np
import pandas as pd
import pytest
from AIS.Channel_Creator import (
    Abstract_Channel_Creator,
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 1024,
}

l_init, l_final, l_step = 350, 1150, 50
wavelength_interval = range(l_init, l_final, l_step)
n = len(wavelength_interval)
ccd_temp = -70


# -----------------------------------------------------------------------------------------------------------
ccd_transmitance_c1 = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 1", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_c1 = np.asarray([float(value) for value in ccd_transmitance_c1])

ccd_transmitance_c2 = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 2", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_c2 = np.asarray([float(value) for value in ccd_transmitance_c2])

ccd_transmitance_c3 = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 3", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_c3 = np.asarray([float(value) for value in ccd_transmitance_c3])

ccd_transmitance_c4 = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 4", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_c4 = np.asarray([float(value) for value in ccd_transmitance_c4])

ccd_transmitance_abs = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 0", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_abs = np.asarray([float(value) for value in ccd_transmitance_abs])


# --------------------------------------------------------------------------------------------------------------------


@pytest.fixture
def abs_chc():
    return Abstract_Channel_Creator(sparc4_operation_mode="phot")


@pytest.fixture
def chc1():
    return Concrete_Channel_1(sparc4_operation_mode="phot")


@pytest.fixture
def chc2():
    return Concrete_Channel_2(sparc4_operation_mode="phot")


@pytest.fixture
def chc3():
    return Concrete_Channel_3(sparc4_operation_mode="phot")


@pytest.fixture
def chc4():
    return Concrete_Channel_4(sparc4_operation_mode="phot")


# ------------------------ Initialize the class --------------------------


def test_sparc_operation_mode(abs_chc):
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
    chc1.calculate_dark_current(ccd_temp)
    assert round(chc1.dark_current, 7) == 5.86e-5


def test_calculate_dark_current_2(chc2):
    chc2.calculate_dark_current(ccd_temp)
    assert chc2.dark_current, 7 == 0.0001467


def test_calculate_dark_current_3(chc3):
    chc3.calculate_dark_current(ccd_temp)
    assert round(chc3.dark_current, 7) == 8.69e-05


def test_calculate_dark_current_4(chc4):
    chc4.calculate_dark_current(ccd_temp)
    assert chc4.dark_current, 7 == 0.0002313


# ----------------------- Calculate read noise  -------------------------


def test_calculate_read_noise_1(chc1):
    chc1.calculate_read_noise(dic)
    assert chc1.read_noise == 6.67


def test_calculate_read_noise_2(chc2):
    chc2.calculate_read_noise(dic)
    assert chc2.read_noise == 6.67


def test_calculate_read_noise_3(chc3):
    chc3.calculate_read_noise(dic)
    assert chc3.read_noise == 6.67


def test_calculate_read_noise_4(chc4):
    chc4.calculate_read_noise(dic)
    assert chc4.read_noise == 6.67


# ------------------- Teste apply_sparc4_spectral_response ----------------------------------------------------------


specific_flux = np.ones((4, n))
specific_flux_c1 = specific_flux[0, :] * 0.5 * ccd_transmitance_c1 / 100
specific_flux_c2 = specific_flux[0, :] * 0.5 * ccd_transmitance_c2 / 100
specific_flux_c3 = specific_flux[0, :] * 0.5 * ccd_transmitance_c3 / 100
specific_flux_c4 = specific_flux[0, :] * 0.5 * ccd_transmitance_c4 / 100

# ---------------------------------------------------------------------------------------------------------------------


# def test_apply_sparc4_spectral_response_photometric_1(chc1):
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc1.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c1)
#     assert np.allclose(specific_extra_ordinary_ray, 0)


# def test_apply_sparc4_spectral_response_polarimetric_1():
#     chc1 = Concrete_Channel_1(sparc4_operation_mode="pol")
#     specific_flux = np.ones((4, n))
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc1.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     specific_flux = specific_flux[0, :] * ccd_transmitance_c1 * 0.5 / 100
#     assert np.allclose(specific_ordinary_ray, specific_flux_c1)
#     assert np.allclose(specific_extra_ordinary_ray, specific_flux_c1)


# def test_apply_sparc4_spectral_response_photometric_2(chc2):
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc2.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c2)
#     assert np.allclose(specific_extra_ordinary_ray, 0)


# def test_apply_sparc4_spectral_response_polarimetric_2():
#     chc2 = Concrete_Channel_2(sparc4_operation_mode="pol")
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc2.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c2)
#     assert np.allclose(specific_extra_ordinary_ray, specific_flux_c2)


# def test_apply_sparc4_spectral_response_photometric_3(chc3):
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc3.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c3)
#     assert np.allclose(specific_extra_ordinary_ray, 0)


# def test_apply_sparc4_spectral_response_polarimetric_3():
#     chc3 = Concrete_Channel_3(sparc4_operation_mode="pol")
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc3.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c3)
#     assert np.allclose(specific_extra_ordinary_ray, specific_flux_c3)


# def test_apply_sparc4_spectral_response_photometric_4(chc4):
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc4.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c4)
#     assert np.allclose(specific_extra_ordinary_ray, 0)


# def test_apply_sparc4_spectral_response_polarimetric_4():
#     chc4 = Concrete_Channel_4(sparc4_operation_mode="pol")
#     (
#         specific_ordinary_ray,
#         specific_extra_ordinary_ray,
#     ) = chc4.apply_sparc4_spectral_response(specific_flux, l_init, l_final, l_step)
#     assert np.allclose(specific_ordinary_ray, specific_flux_c4)
#     assert np.allclose(specific_extra_ordinary_ray, specific_flux_c4)
