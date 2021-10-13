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
from Spectrum_Calculation import Spectrum_Calculation

from tests.test_AIS_operation import multiply_matrices

from .SPARC4_SR_curves import (
    analyser_extra_ordinary_ray,
    analyser_ordinary_ray,
    calibration_wheel,
    camera_c1,
    camera_c2,
    camera_c3,
    camera_c4,
    ccd_transmitance_c1,
    ccd_transmitance_c2,
    ccd_transmitance_c3,
    ccd_transmitance_c4,
    colimator_transmitance,
    dichroic_c1,
    dichroic_c2,
    dichroic_c3,
    dichroic_c4,
    retarder,
    specific_flux,
    wavelength_interval,
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

magnitude = 22
star_temperature = 5700
ccd_temp = -70


# --------------------------------------------------------------------------------------------------------------------


@pytest.fixture
def abs_chc():
    return Abstract_Channel_Creator(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


@pytest.fixture
def chc1():
    return Concrete_Channel_1(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


@pytest.fixture
def chc2():
    return Concrete_Channel_2(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


@pytest.fixture
def chc3():
    return Concrete_Channel_3(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


@pytest.fixture
def chc4():
    return Concrete_Channel_4(
        sparc4_operation_mode="phot", wavelength_interval=wavelength_interval
    )


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


# ------------------- Teste apply_sparc4_spectral_response -----------------------------

specific_flux_c1 = multiply_matrices(calibration_wheel, specific_flux.copy())
specific_flux_c1 = multiply_matrices(retarder, specific_flux_c1)
c1_ordinary_ray = multiply_matrices(analyser_ordinary_ray, specific_flux_c1.copy())
c1_extra_ordinary_ray = multiply_matrices(analyser_extra_ordinary_ray, specific_flux_c1)
c1_ordinary_ray = np.multiply(colimator_transmitance, c1_ordinary_ray)
c1_extra_ordinary_ray = np.multiply(colimator_transmitance, c1_extra_ordinary_ray)
c1_ordinary_ray = np.multiply(dichroic_c1, c1_ordinary_ray)
c1_extra_ordinary_ray = np.multiply(dichroic_c1, c1_extra_ordinary_ray)
c1_ordinary_ray = np.multiply(camera_c1, c1_ordinary_ray)
c1_extra_ordinary_ray = np.multiply(camera_c1, c1_extra_ordinary_ray)
c1_ordinary_ray = np.multiply(ccd_transmitance_c1, c1_ordinary_ray)
c1_extra_ordinary_ray = np.multiply(ccd_transmitance_c1, c1_extra_ordinary_ray)


specific_flux_c2 = multiply_matrices(calibration_wheel, specific_flux.copy())
specific_flux_c2 = multiply_matrices(retarder, specific_flux_c2)
c2_ordinary_ray = multiply_matrices(analyser_ordinary_ray, specific_flux_c2.copy())
c2_extra_ordinary_ray = multiply_matrices(analyser_extra_ordinary_ray, specific_flux_c2)
c2_ordinary_ray = np.multiply(colimator_transmitance, c2_ordinary_ray)
c2_extra_ordinary_ray = np.multiply(colimator_transmitance, c2_extra_ordinary_ray)
c2_ordinary_ray = np.multiply(dichroic_c2, c2_ordinary_ray)
c2_extra_ordinary_ray = np.multiply(dichroic_c2, c2_extra_ordinary_ray)
c2_ordinary_ray = np.multiply(camera_c2, c2_ordinary_ray)
c2_extra_ordinary_ray = np.multiply(camera_c2, c2_extra_ordinary_ray)
c2_ordinary_ray = np.multiply(ccd_transmitance_c2, c2_ordinary_ray)
c2_extra_ordinary_ray = np.multiply(ccd_transmitance_c2, c2_extra_ordinary_ray)

specific_flux_c3 = multiply_matrices(calibration_wheel, specific_flux.copy())
specific_flux_c3 = multiply_matrices(retarder, specific_flux_c3)
c3_ordinary_ray = multiply_matrices(analyser_ordinary_ray, specific_flux_c3.copy())
c3_extra_ordinary_ray = multiply_matrices(analyser_extra_ordinary_ray, specific_flux_c3)
c3_ordinary_ray = np.multiply(colimator_transmitance, c3_ordinary_ray)
c3_extra_ordinary_ray = np.multiply(colimator_transmitance, c3_extra_ordinary_ray)
c3_ordinary_ray = np.multiply(dichroic_c3, c3_ordinary_ray)
c3_extra_ordinary_ray = np.multiply(dichroic_c3, c3_extra_ordinary_ray)
c3_ordinary_ray = np.multiply(camera_c3, c3_ordinary_ray)
c3_extra_ordinary_ray = np.multiply(camera_c3, c3_extra_ordinary_ray)
c3_ordinary_ray = np.multiply(ccd_transmitance_c3, c3_ordinary_ray)
c3_extra_ordinary_ray = np.multiply(ccd_transmitance_c3, c3_extra_ordinary_ray)

specific_flux_c4 = multiply_matrices(calibration_wheel, specific_flux.copy())
specific_flux_c4 = multiply_matrices(retarder, specific_flux_c4)
c4_ordinary_ray = multiply_matrices(analyser_ordinary_ray, specific_flux_c4.copy())
c4_extra_ordinary_ray = multiply_matrices(analyser_extra_ordinary_ray, specific_flux_c4)
c4_ordinary_ray = np.multiply(colimator_transmitance, c4_ordinary_ray)
c4_extra_ordinary_ray = np.multiply(colimator_transmitance, c4_extra_ordinary_ray)
c4_ordinary_ray = np.multiply(dichroic_c4, c4_ordinary_ray)
c4_extra_ordinary_ray = np.multiply(dichroic_c4, c4_extra_ordinary_ray)
c4_ordinary_ray = np.multiply(camera_c4, c4_ordinary_ray)
c4_extra_ordinary_ray = np.multiply(camera_c4, c4_extra_ordinary_ray)
c4_ordinary_ray = np.multiply(ccd_transmitance_c4, c4_ordinary_ray)
c4_extra_ordinary_ray = np.multiply(ccd_transmitance_c4, c4_extra_ordinary_ray)

# ---------------------------------------------------------------------------------------------------------------------


def test_apply_sparc4_spectral_response_polarimetric_1():
    chc1_pol = Concrete_Channel_1(
        sparc4_operation_mode="pol", wavelength_interval=wavelength_interval
    )
    (
        specific_star_ordinary_ray,
        specific_star_extra_ordinary_ray,
    ) = chc1_pol.apply_sparc4_spectral_response(specific_flux.copy())
    assert np.allclose(specific_star_ordinary_ray, c1_ordinary_ray)
    assert np.allclose(specific_star_extra_ordinary_ray, c1_extra_ordinary_ray)


def test_apply_sparc4_spectral_response_polarimetric_2():
    chc2_pol = Concrete_Channel_2(
        sparc4_operation_mode="pol", wavelength_interval=wavelength_interval
    )
    (
        specific_star_ordinary_ray,
        specific_star_extra_ordinary_ray,
    ) = chc2_pol.apply_sparc4_spectral_response(
        specific_flux.copy(),
    )
    assert np.allclose(specific_star_ordinary_ray, c2_ordinary_ray)
    assert np.allclose(specific_star_extra_ordinary_ray, c2_extra_ordinary_ray)


def test_apply_sparc4_spectral_response_polarimetric_3():
    chc3_pol = Concrete_Channel_3(
        sparc4_operation_mode="pol", wavelength_interval=wavelength_interval
    )
    (
        specific_star_ordinary_ray,
        specific_star_extra_ordinary_ray,
    ) = chc3_pol.apply_sparc4_spectral_response(
        specific_flux.copy(),
    )
    assert np.allclose(specific_star_ordinary_ray, c3_ordinary_ray)
    assert np.allclose(specific_star_extra_ordinary_ray, c3_extra_ordinary_ray)


def test_apply_sparc4_spectral_response_polarimetric_4():
    chc4_pol = Concrete_Channel_4(
        sparc4_operation_mode="pol", wavelength_interval=wavelength_interval
    )
    (
        specific_star_ordinary_ray,
        specific_star_extra_ordinary_ray,
    ) = chc4_pol.apply_sparc4_spectral_response(
        specific_flux.copy(),
    )
    assert np.allclose(specific_star_ordinary_ray, c4_ordinary_ray)
    assert np.allclose(specific_star_extra_ordinary_ray, c4_extra_ordinary_ray)


# -------------------------------------------------------------------------------------------------------------------------
c1_ordinary_ray_phot = np.multiply(colimator_transmitance, specific_flux.copy())
c1_ordinary_ray_phot = np.multiply(dichroic_c1, c1_ordinary_ray_phot)
c1_ordinary_ray_phot = np.multiply(camera_c1, c1_ordinary_ray_phot)
c1_ordinary_ray_phot = np.multiply(ccd_transmitance_c1, c1_ordinary_ray_phot)
c1_extra_ordinary_ray_phot = 0

c2_ordinary_ray_phot = np.multiply(colimator_transmitance, specific_flux.copy())
c2_ordinary_ray_phot = np.multiply(dichroic_c2, c2_ordinary_ray_phot)
c2_ordinary_ray_phot = np.multiply(camera_c2, c2_ordinary_ray_phot)
c2_ordinary_ray_phot = np.multiply(ccd_transmitance_c2, c2_ordinary_ray_phot)
c2_extra_ordinary_ray_phot = 0

c3_ordinary_ray_phot = np.multiply(colimator_transmitance, specific_flux.copy())
c3_ordinary_ray_phot = np.multiply(dichroic_c3, c3_ordinary_ray_phot)
c3_ordinary_ray_phot = np.multiply(camera_c3, c3_ordinary_ray_phot)
c3_ordinary_ray_phot = np.multiply(ccd_transmitance_c3, c3_ordinary_ray_phot)
c3_extra_ordinary_ray_phot = 0

c4_ordinary_ray_phot = np.multiply(colimator_transmitance, specific_flux.copy())
c4_ordinary_ray_phot = np.multiply(dichroic_c4, c4_ordinary_ray_phot)
c4_ordinary_ray_phot = np.multiply(camera_c4, c4_ordinary_ray_phot)
c4_ordinary_ray_phot = np.multiply(ccd_transmitance_c4, c4_ordinary_ray_phot)
c4_extra_ordinary_ray_phot = 0


def test_apply_sparc4_spectral_response_photometric_1(chc1):
    star_specific_photons_per_second = chc1.apply_sparc4_spectral_response(
        specific_flux.copy(),
    )
    assert np.allclose(star_specific_photons_per_second, c1_ordinary_ray_phot)


def test_apply_sparc4_spectral_response_photometric_2(chc2):
    star_specific_photons_per_second = chc2.apply_sparc4_spectral_response(
        specific_flux.copy()
    )
    assert np.allclose(star_specific_photons_per_second, c2_ordinary_ray_phot)


def test_apply_sparc4_spectral_response_photometric_3(chc3):

    star_specific_photons_per_second = chc3.apply_sparc4_spectral_response(
        specific_flux.copy()
    )
    assert np.allclose(star_specific_photons_per_second, c3_ordinary_ray_phot)


def test_apply_sparc4_spectral_response_photometric_4(chc4):

    star_specific_photons_per_second = chc4.apply_sparc4_spectral_response(
        specific_flux.copy()
    )
    assert np.allclose(star_specific_photons_per_second, c4_ordinary_ray_phot)
