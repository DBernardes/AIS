# -*- coding: utf-8 -*-
"""AIS operation tests.

This script presents the tests of the AIS operation. These tests are:

    - Create Image Name: tess of the function that creates the image name based
    on the provided operation mode of the CCD

    - Configura Gain: tess of the function that sets the CCD gain based
    on the provided operation mode of the CCD


Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


import os

import numpy as np
import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

from .AIS_spectral_response_curves import (
    POLARIZER_90_MATRIX,
    POLARIZER_MATRIX,
    THETA_POL,
    air_mass,
    analyser_transmitance,
    atm_transmitance,
    calculate_retarder_matrix,
    camera_c1,
    ccd_operation_mode,
    ccd_transmitance_c1,
    collimator_transmitance,
    l_final,
    l_init,
    l_step,
    multiply_matrices,
    polarizer_transmitance,
    retardance_quarter,
    retarder_transmitance,
    sky_condition,
    sky_specific_photons_per_second,
    sparc4_operation_mode,
    star_specific_photons_per_second,
    star_temperature,
    tel_reflectance,
)

star_specific_photons_per_second = star_specific_photons_per_second.copy()
channel = 1
star_coordinates = (100, 100)
image_dir = "a"
bias_level = 500


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=ccd_operation_mode,
        channel=channel,
        star_coordinates=star_coordinates,
        bias_level=bias_level,
        sparc4_operation_mode=sparc4_operation_mode,
        image_dir=image_dir,
        wavelength_interval=(l_init, l_final, l_step),
        star_temperature=star_temperature,
        air_mass=air_mass,
        sky_condition=sky_condition,
    )


# ----------------------- Calculate dark current  -------------------------


def test_calculate_dark_current(ais):
    ais._calculate_dark_current()
    assert round(ais.dark_current, 7) == 5.86e-5


# -------------------------Calculate Read Noise -------------------------


def test_calculate_read_noise(ais):
    ais._calculate_read_noise(ccd_operation_mode)
    assert ais.read_noise == 6.67


# ----------------------- Apply spectruns ---------------------------


def test_calculate_star_specific_photons_per_second(ais):
    ais._calculate_star_specific_photons_per_second()
    assert np.allclose(
        ais.star_specific_photons_per_second, star_specific_photons_per_second.copy()
    )


def test_calculate_sky_specific_photons_per_second(ais):
    ais._calculate_sky_specific_photons_per_second()
    assert np.allclose(
        ais.sky_specific_photons_per_second, sky_specific_photons_per_second.copy()
    )


def test_apply_atmosphere_spectral_response_star(ais):
    new_star_specific_photons_per_second = np.multiply(
        star_specific_photons_per_second[0][0].copy(), atm_transmitance
    )
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(
        ais.star_specific_photons_per_second[0][0][0],
        new_star_specific_photons_per_second,
    )


def test_apply_atmosphere_spectral_response_sky(ais):
    new_sky_specific_photons_per_second = np.multiply(
        sky_specific_photons_per_second[0][0].copy(), atm_transmitance
    )
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(
        ais.sky_specific_photons_per_second[0][0][0].copy(),
        new_sky_specific_photons_per_second,
    )


def test_apply_telescope_spectral_response_star(ais):
    ais.apply_telescope_spectral_response()
    new_star_specific_photons_per_second = np.multiply(
        star_specific_photons_per_second.copy(), tel_reflectance
    )
    assert np.allclose(
        ais.star_specific_photons_per_second, new_star_specific_photons_per_second
    )


def test_apply_telescope_spectral_response_sky(ais):
    new_sky_specific_photons_per_second = np.multiply(
        sky_specific_photons_per_second.copy(), tel_reflectance
    )
    ais.apply_telescope_spectral_response()
    assert np.allclose(
        ais.sky_specific_photons_per_second, new_sky_specific_photons_per_second
    )


# -----------------------------------Test apply SPARC4 spectral response ----------------------------------------------------


# def test_apply_sparc4_spectral_response_star():
#     star_specific_photons_per_second_c1 = star_specific_photons_per_second.copy()
#     # Applying collimator
#     star_specific_photons_per_second_c1[0][0] = np.multiply(
#         collimator_transmitance, star_specific_photons_per_second_c1[0][0]
#     )
#     # Applying camera
#     star_specific_photons_per_second_c1[0][0] = np.multiply(
#         camera_c1, star_specific_photons_per_second_c1[0][0]
#     )
#     # Applying ccd
#     star_specific_photons_per_second_c1[0][0] = np.multiply(
#         ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
#     )

#     sparc4_operation_mode = {
#         "acquisition_mode": "photometric",
#     }
#     ais = Artificial_Image_Simulator(
#         ccd_operation_mode=ccd_operation_mode,
#         sparc4_operation_mode=sparc4_operation_mode,
#         channel=1,
#     )
#     ais.apply_sparc4_spectral_response()
#     assert np.allclose(
#         star_specific_photons_per_second_c1, ais.star_specific_photons_per_second
#     )


# ----------------------------------------------------------------------------------------------------------------


# def test_apply_sparc4_spectral_response_sky(ais):
#     ais.apply_sparc4_spectral_response()

#     assert np.allclose(ais.sky_specific_photons_per_second, new_ordinary_ray)


# -----------------------------test _integrate_fluxes ------------------------


def test_integrate_specific_photons_per_second(ais):
    photons_per_second = np.trapz(star_specific_photons_per_second[0][0, :])
    ais.star_specific_photons_per_second = star_specific_photons_per_second
    ais.sky_specific_photons_per_second = star_specific_photons_per_second
    ais._integrate_specific_photons_per_second()
    assert ais.star_ordinary_ray == photons_per_second
    assert ais.sky_photons_per_second == photons_per_second


# -----------------------------test _create_image_name------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, t_exp, image_name",
    [
        ("Conv", 1, 0.1, 1, 1, 1, "CONV_HSS0.1_PA1_B1_TEXP1_G1"),
        ("Conv", 1, 0.1, 1, 2, 1, "CONV_HSS0.1_PA1_B2_TEXP1_G1"),
        ("Conv", 1, 0.1, 2, 1, 1, "CONV_HSS0.1_PA2_B1_TEXP1_G1"),
        ("Conv", 1, 0.1, 2, 2, 1, "CONV_HSS0.1_PA2_B2_TEXP1_G1"),
        ("Conv", 1, 1, 1, 1, 1, "CONV_HSS1_PA1_B1_TEXP1_G1"),
        ("Conv", 1, 1, 1, 2, 1, "CONV_HSS1_PA1_B2_TEXP1_G1"),
        ("Conv", 1, 1, 2, 1, 1, "CONV_HSS1_PA2_B1_TEXP1_G1"),
        ("Conv", 1, 1, 2, 2, 1, "CONV_HSS1_PA2_B2_TEXP1_G1"),
        ("EM", 2, 1, 1, 1, 1, "EM_HSS1_PA1_B1_TEXP1_G2"),
        ("EM", 2, 1, 1, 2, 1, "EM_HSS1_PA1_B2_TEXP1_G2"),
        ("EM", 2, 1, 2, 1, 1, "EM_HSS1_PA2_B1_TEXP1_G2"),
        ("EM", 2, 1, 2, 2, 1, "EM_HSS1_PA2_B2_TEXP1_G2"),
        ("EM", 2, 10, 1, 1, 1, "EM_HSS10_PA1_B1_TEXP1_G2"),
        ("EM", 2, 10, 1, 2, 1, "EM_HSS10_PA1_B2_TEXP1_G2"),
        ("EM", 2, 10, 2, 1, 1, "EM_HSS10_PA2_B1_TEXP1_G2"),
        ("EM", 2, 10, 2, 2, 1, "EM_HSS10_PA2_B2_TEXP1_G2"),
        ("EM", 2, 20, 1, 1, 1, "EM_HSS20_PA1_B1_TEXP1_G2"),
        ("EM", 2, 20, 1, 2, 1, "EM_HSS20_PA1_B2_TEXP1_G2"),
        ("EM", 2, 20, 2, 1, 1, "EM_HSS20_PA2_B1_TEXP1_G2"),
        ("EM", 2, 20, 2, 2, 1, "EM_HSS20_PA2_B2_TEXP1_G2"),
        ("EM", 2, 30, 1, 1, 1, "EM_HSS30_PA1_B1_TEXP1_G2"),
        ("EM", 2, 30, 1, 2, 1, "EM_HSS30_PA1_B2_TEXP1_G2"),
        ("EM", 2, 30, 2, 1, 1, "EM_HSS30_PA2_B1_TEXP1_G2"),
        ("EM", 2, 30, 2, 2, 1, "EM_HSS30_PA2_B2_TEXP1_G2"),
        ("EM", 2, 30, 2, 2, 2, "EM_HSS30_PA2_B2_TEXP2_G2"),
    ],
)
def test_create_image_name(ais, em_mode, em_gain, hss, preamp, binn, t_exp, image_name):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": t_exp,
        "ccd_temp": -70,
        "image_size": 200,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode)
    ais._configure_image_name()
    assert ais.image_name == image_name


# -----------------------------test _configure_gain------------------------


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
def test_configure_gain(ais, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    ccd_operation_mode = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    ais.ccd_operation_mode = ccd_operation_mode
    ais._configure_gain()
    assert ais.ccd_gain == ccd_gain


# --------------------------- test create artificial image ----------------


def test_create_artificial_image_phot():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_artificial_image()


def test_create_artificial_image_pol():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode,
        image_dir=os.path.join("FITS"),
        sparc4_operation_mode=sparc4_operation_mode,
    )
    ais.apply_sparc4_spectral_response()
    ais.create_artificial_image()


def test_create_background_image():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_background_image()


def test_creat_bias_image():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_bias_image()


def test_creat_dark_image():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_dark_image()


def test_creat_random_image():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_random_image(n=2)


def test_creat_flat_image():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_flat_image()