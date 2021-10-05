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
from scipy.interpolate import splev, splrep
from Spectrum_Calculation import Spectrum_Calculation

from .SPARC4_SR_curves import (
    analyser_extra_ordinary_ray,
    analyser_ordinary_ray,
    calibration_wheel,
    camera_c1,
    ccd_transmitance_c1,
    colimator_transmitance,
    dichroic_c1,
    l_final,
    l_init,
    l_step,
    retarder,
    wavelength_interval,
    wavelength_interval_len,
)

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 100,
}

channel = 1
star_coordinates = (100, 100)
sparc4_operation_mode = "pol"
image_dir = "a"
star_temperature = 5700
star_magnitude = 22
bias_level = 500
sc = Spectrum_Calculation(wavelength_interval, star_temperature)
star_specific_flux = sc.calculate_specific_flux(star_magnitude)
sky_specific_flux = sc.calculate_specific_flux(star_magnitude + 3)
specific_flux_length = len(star_specific_flux)


def multiply_matrices(matrix, specific_flux):
    for i in range(len(specific_flux[0])):
        specific_flux[:, i] = np.dot(matrix, specific_flux[:, i])
    return specific_flux


# ----------------------- importing the telescope spectral response ----------------
ss = pd.read_csv(
    os.path.join("Telescope_Spectral_Response", "telescope_spectral_response.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

tel_wavelength_interval = ss["(nm)"]
tel_reflectance = ss["(%)"] / 100
spl = splrep(tel_wavelength_interval, tel_reflectance)
tel_reflectance = splev(wavelength_interval, spl)
# ------------------------------------------------------------------------------------


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=dic,
        channel=channel,
        star_coordinates=star_coordinates,
        bias_level=bias_level,
        sparc4_operation_mode=sparc4_operation_mode,
        image_dir=image_dir,
        wavelength_interval=(l_init, l_final, l_step),
        star_temperature=star_temperature,
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
    assert np.allclose(ais.specific_star_ordinary_ray, star_specific_flux)
    assert np.allclose(
        ais.specific_star_extra_ordinary_ray, np.zeros((4, wavelength_interval_len))
    )


def test_calculate_sky_specific_flux(ais):
    ais._calculate_sky_specific_flux()
    assert np.allclose(ais.specific_sky_ordinary_ray, sky_specific_flux.copy())
    assert np.allclose(
        ais.specific_sky_extra_ordinary_ray, np.zeros((4, wavelength_interval_len))
    )


def test_apply_atmosphere_spectral_response_star(ais):
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(ais.specific_star_ordinary_ray, star_specific_flux)


def test_apply_atmosphere_spectral_response_sky(ais):
    ais.apply_atmosphere_spectral_response()
    assert np.allclose(ais.specific_sky_ordinary_ray, sky_specific_flux.copy())


def test_apply_telescope_spectral_response_star(ais):
    ais.apply_telescope_spectral_response()
    new_star_specific_flux = np.multiply(star_specific_flux, tel_reflectance)
    assert np.allclose(ais.specific_star_ordinary_ray, new_star_specific_flux)


def test_apply_telescope_spectral_response_sky(ais):
    new_sky_specific_flux = np.multiply(sky_specific_flux.copy(), tel_reflectance)
    ais.apply_telescope_spectral_response()
    assert np.allclose(ais.specific_sky_ordinary_ray, new_sky_specific_flux)


def test_apply_sparc4_spectral_response_star(ais):

    ais.apply_sparc4_spectral_response()
    specific_flux = multiply_matrices(calibration_wheel, star_specific_flux.copy())
    specific_flux = multiply_matrices(retarder, specific_flux)
    new_ordinary_ray = multiply_matrices(analyser_ordinary_ray, specific_flux.copy())
    new_extra_ordinary_ray = multiply_matrices(
        analyser_extra_ordinary_ray, specific_flux.copy()
    )
    new_ordinary_ray = np.multiply(colimator_transmitance, new_ordinary_ray)
    new_extra_ordinary_ray = np.multiply(colimator_transmitance, new_extra_ordinary_ray)
    new_ordinary_ray = np.multiply(new_ordinary_ray, dichroic_c1)
    new_extra_ordinary_ray = np.multiply(new_extra_ordinary_ray, dichroic_c1)
    new_ordinary_ray = np.multiply(new_ordinary_ray, ccd_transmitance_c1)
    new_extra_ordinary_ray = np.multiply(new_extra_ordinary_ray, ccd_transmitance_c1)
    assert np.allclose(ais.specific_star_ordinary_ray, new_ordinary_ray)
    assert np.allclose(ais.specific_star_extra_ordinary_ray, new_extra_ordinary_ray)


def test_apply_sparc4_spectral_response_sky(ais):
    ais.apply_sparc4_spectral_response()

    new_sky_specific_flux = multiply_matrices(calibration_wheel, sky_specific_flux)
    new_sky_specific_flux = multiply_matrices(retarder, new_sky_specific_flux)
    new_ordinary_ray = multiply_matrices(
        analyser_ordinary_ray, new_sky_specific_flux.copy()
    )

    new_extra_ordinary_ray = multiply_matrices(
        analyser_extra_ordinary_ray, new_sky_specific_flux
    )
    new_ordinary_ray = np.multiply(colimator_transmitance, new_ordinary_ray)
    new_extra_ordinary_ray = np.multiply(colimator_transmitance, new_extra_ordinary_ray)
    new_ordinary_ray = np.multiply(new_ordinary_ray, dichroic_c1)
    new_extra_ordinary_ray = np.multiply(new_extra_ordinary_ray, dichroic_c1)
    new_ordinary_ray = np.multiply(new_ordinary_ray, camera_c1)
    new_extra_ordinary_ray = np.multiply(new_extra_ordinary_ray, camera_c1)
    new_ordinary_ray = np.multiply(new_ordinary_ray, ccd_transmitance_c1)
    new_extra_ordinary_ray = np.multiply(new_extra_ordinary_ray, ccd_transmitance_c1)

    assert np.allclose(ais.specific_sky_ordinary_ray, new_ordinary_ray)
    assert np.allclose(ais.specific_sky_extra_ordinary_ray, new_extra_ordinary_ray)


# -----------------------------test _integrate_fluxes ------------------------


def test_integrate_fluxes(ais):
    photons_per_second = np.trapz(star_specific_flux[0, :])
    ais.specific_star_ordinary_ray = star_specific_flux
    ais.specific_star_extra_ordinary_ray = np.zeros((4, wavelength_interval_len))
    ais.specific_sky_ordinary_ray = star_specific_flux
    ais.specific_sky_extra_ordinary_ray = np.zeros((4, wavelength_interval_len))
    ais._integrate_specific_fluxes()
    assert ais.star_ordinary_ray == photons_per_second
    assert ais.star_extra_ordinary_ray == 0
    assert ais.sky_ordinary_ray == photons_per_second
    assert ais.sky_extra_ordinary_ray == 0


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
        "image_size": 100,
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
        "image_size": 100,
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
        "image_size": 100,
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
        "image_size": 100,
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
        "image_size": 100,
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
        "image_size": 100,
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
        "image_size": 100,
    }
    ais = Artificial_Image_Simulator(dic, image_dir=os.path.join("..", "FITS"))
    ais.apply_sparc4_spectral_response()
    ais.create_flat_image()
