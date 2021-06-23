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
import pandas as pd
import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from Spectrum_Calculation import Spectrum_Calculation

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 200,
}

sc = Spectrum_Calculation(5700, 350, 1150, 50)
star_specific_flux = sc.calculate_star_specific_flux()
sky_specific_flux = sc.calculate_sky_specific_flux()
specific_flux_length = len(star_specific_flux)
wavelength_interval = range(350, 1150, 50)
ccd_transmitance_c1 = np.asarray(
    pd.read_excel(os.path.join("SPARC4_Spectral_Response", "Channel 1", "ccd.xlsx"))
)[1:, 1]
ccd_transmitance_c1 = np.asarray([float(value) for value in ccd_transmitance_c1])


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=dic,
        channel=1,
        gaussian_std=3,
        star_coordinates=(100, 100),
        bias_level=500,
        sparc4_operation_mode="pol",
        image_dir="a",
        star_wavelength_interval=(350, 1150, 50),
        star_temperature=5700,
    )


# ----------------------- Calculate dark current  -------------------------


def test_calculate_dark_current(ais):
    ais._calculate_dark_current()
    assert round(ais.dark_current, 7) == 5.86e-5


# -------------------------Calculate Read Noise -------------------------


def test_calculate_read_noise(ais):
    ais._calculate_read_noise(dic)
    assert ais.read_noise == 6.67


# ----------------------- Apply spectruns ---------------------------


def test_calculate_star_specific_flux(ais):
    ais._calculate_star_specific_flux()
    boolean_test = ais.star_specific_flux == star_specific_flux
    assert boolean_test.all()


def test_calculate_sky_specific_flux(ais):
    ais._calculate_sky_specific_flux()
    boolean_test = ais.sky_specific_flux == sky_specific_flux
    assert boolean_test.all()


def test_apply_atmosphere_spectral_response_star(ais):
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(ais.star_specific_flux, star_specific_flux)


# def test_apply_atmosphere_spectral_response_sky(ais):
#     ais.apply_atmosphere_spectral_response()
#     assert np.allclose(ais.sky_specific_flux, sky_specific_flux)


def test_apply_telescope_spectral_response_star(ais):
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(ais.star_specific_flux, star_specific_flux)


# def test_apply_telescope_spectral_response_sky(ais):
#     ais.apply_atmosphere_spectral_response()
#     assert np.allclose(ais.sky_specific_flux, sky_specific_flux)


# def test_apply_sparc4_spectral_response_star(ais):
#     specific_flux = (
#         np.multiply(star_specific_flux[0, :], ccd_transmitance_c1) * 0.5 / 100
#     )
#     ais.apply_sparc4_spectral_response()
#     assert np.allclose(ais.specific_ordinary_ray, specific_flux)
#     assert np.allclose(ais.specific_extra_ordinary_ray, specific_flux)


# def test_apply_sparc4_spectral_response_sky(ais):
#     specific_flux = np.multiply(sky_specific_flux[0, :], ccd_transmitance_c1) / 100
#     ais.apply_sparc4_spectral_response()
#     assert np.allclose(ais.ordinary_ray, specific_flux)


# -----------------------------test _create_image_name------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, t_exp, image_name",
    [
        (0, 1, 0.1, 1, 1, 1, "CONV_HSS0.1_PA1_B1_TEXP1_G1"),
        (0, 1, 0.1, 1, 2, 1, "CONV_HSS0.1_PA1_B2_TEXP1_G1"),
        (0, 1, 0.1, 2, 1, 1, "CONV_HSS0.1_PA2_B1_TEXP1_G1"),
        (0, 1, 0.1, 2, 2, 1, "CONV_HSS0.1_PA2_B2_TEXP1_G1"),
        (0, 1, 1, 1, 1, 1, "CONV_HSS1_PA1_B1_TEXP1_G1"),
        (0, 1, 1, 1, 2, 1, "CONV_HSS1_PA1_B2_TEXP1_G1"),
        (0, 1, 1, 2, 1, 1, "CONV_HSS1_PA2_B1_TEXP1_G1"),
        (0, 1, 1, 2, 2, 1, "CONV_HSS1_PA2_B2_TEXP1_G1"),
        (1, 2, 1, 1, 1, 1, "EM_HSS1_PA1_B1_TEXP1_G2"),
        (1, 2, 1, 1, 2, 1, "EM_HSS1_PA1_B2_TEXP1_G2"),
        (1, 2, 1, 2, 1, 1, "EM_HSS1_PA2_B1_TEXP1_G2"),
        (1, 2, 1, 2, 2, 1, "EM_HSS1_PA2_B2_TEXP1_G2"),
        (1, 2, 10, 1, 1, 1, "EM_HSS10_PA1_B1_TEXP1_G2"),
        (1, 2, 10, 1, 2, 1, "EM_HSS10_PA1_B2_TEXP1_G2"),
        (1, 2, 10, 2, 1, 1, "EM_HSS10_PA2_B1_TEXP1_G2"),
        (1, 2, 10, 2, 2, 1, "EM_HSS10_PA2_B2_TEXP1_G2"),
        (1, 2, 20, 1, 1, 1, "EM_HSS20_PA1_B1_TEXP1_G2"),
        (1, 2, 20, 1, 2, 1, "EM_HSS20_PA1_B2_TEXP1_G2"),
        (1, 2, 20, 2, 1, 1, "EM_HSS20_PA2_B1_TEXP1_G2"),
        (1, 2, 20, 2, 2, 1, "EM_HSS20_PA2_B2_TEXP1_G2"),
        (1, 2, 30, 1, 1, 1, "EM_HSS30_PA1_B1_TEXP1_G2"),
        (1, 2, 30, 1, 2, 1, "EM_HSS30_PA1_B2_TEXP1_G2"),
        (1, 2, 30, 2, 1, 1, "EM_HSS30_PA2_B1_TEXP1_G2"),
        (1, 2, 30, 2, 2, 1, "EM_HSS30_PA2_B2_TEXP1_G2"),
        (1, 2, 30, 2, 2, 2, "EM_HSS30_PA2_B2_TEXP2_G2"),
    ],
)
def test_create_image_name(ais, em_mode, em_gain, hss, preamp, binn, t_exp, image_name):
    dic = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": t_exp,
        "ccd_temp": -70,
        "image_size": 200,
    }
    ais = Artificial_Image_Simulator(dic)
    ais._configure_image_name()
    assert ais.image_name == image_name


# -----------------------------test _configure_gain------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, ccd_gain",
    [
        (0, 2, 0.1, 1, 1, 3.35),
        (0, 2, 0.1, 2, 1, 0.80),
        (0, 2, 1, 1, 1, 3.37),
        (0, 2, 1, 2, 1, 0.80),
        (1, 2, 1, 1, 1, 15.90),
        (1, 2, 1, 2, 1, 3.88),
        (1, 2, 10, 1, 1, 16.00),
        (1, 2, 10, 2, 1, 3.96),
        (1, 2, 20, 1, 1, 16.40),
        (1, 2, 20, 2, 1, 4.39),
        (1, 2, 30, 1, 1, 17.20),
        (1, 2, 30, 2, 1, 5.27),
    ],
)
def test_configure_gain(ais, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    dic = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    ais.ccd_operation_mode = dic
    ais._configure_gain()
    assert ais.ccd_gain == ccd_gain


# --------------------------- test create artificial image ----------------


def test_create_artificial_image_phot():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_artificial_image()


def test_create_artificial_image_pol():
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
    ais = Artificial_Image_Simulator(
        dic, image_dir=os.path.join("..", "FITS"), sparc4_operation_mode="pol"
    )
    ais.apply_sparc4_spectral_response()
    ais.create_artificial_image()


def test_create_background_image():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_background_image()


def test_creat_bias_image():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_bias_image()


def test_creat_dark_image():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_dark_image()


def test_creat_random_image():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_random_image(n=2)


def test_creat_flat_image():
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
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_flat_image()
