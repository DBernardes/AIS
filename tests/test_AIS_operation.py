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
import datetime
import numpy as np
import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Spectral_Energy_Distribution import Source, Sky
from AIS.Spectral_Response import Atmosphere, Telescope, Channel

from .AIS_spectral_response_curves import ccd_operation_mode

calculation_method = 'blackbody'
magnitude = 10
wavelegnth_interval = (350, 1100, 100)
temperature = 5700
channel_id = 1


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode, channel_id, temperature)


# ------------------------------------------------------------


def test_print_available_spectral_types(ais):
    ais.print_available_spectral_types()


def test_create_source_sed_blackbody(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    src = Source()
    wv2, sed2 = src.calculate_sed(
        calculation_method, magnitude, wavelegnth_interval, temperature)
    assert np.allclose(ais.wavelength, wv2)
    assert np.allclose(ais.source_sed, sed2)


def test_create_source_sed_spectral_lib(ais):
    calculation_method = 'spectral_library'
    spectral_type = 'A0V'
    ais.create_source_sed(
        calculation_method, magnitude, spectral_type=spectral_type)
    src = Source()
    wv2, sed2 = src.calculate_sed(
        calculation_method, magnitude, spectral_type=spectral_type)
    assert np.allclose(ais.wavelength, wv2)
    assert np.allclose(ais.source_sed, sed2)


wv = wavelegnth_interval
obj_wavelength = np.linspace(wv[0], wv[1], wv[2])
moon_phase = 'new'


def test_create_sky_sed(ais):
    ais.create_source_sed('blackbody', magnitude,
                          wavelegnth_interval, temperature)
    ais.create_sky_sed(moon_phase)
    sky = Sky()
    sed2 = sky.calculate_sed(moon_phase, obj_wavelength)

    assert np.allclose(ais.sky_sed, sed2)


# # ----------------------- Apply spectruns ---------------------------
air_mass = 1
sky_condition = 'photometric'


def test_apply_atmosphere_spectral_response(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    atm = Atmosphere()
    new_sed = atm.apply_spectral_response(
        ais.source_sed, obj_wavelength, air_mass, sky_condition)
    ais.apply_atmosphere_spectral_response(air_mass, sky_condition)
    assert np.allclose(ais.source_sed, new_sed)


def test_apply_telescope_spectral_response(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    ais.create_sky_sed(moon_phase)
    tel = Telescope()
    new_sed = tel.apply_spectral_response(ais.source_sed, obj_wavelength)
    new_sky_sed = tel.apply_spectral_response(ais.sky_sed, obj_wavelength)
    ais.apply_telescope_spectral_response()
    assert np.allclose(ais.source_sed, new_sed)
    assert np.allclose(ais.sky_sed, new_sky_sed)


# -----------------------------------Test apply SPARC4 spectral response ----------------------------------------------------


def test_apply_sparc4_spectral_response_photometric(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    ais.create_sky_sed(moon_phase)
    channel = Channel(channel_id)
    channel.write_sparc4_operation_mode('photometric')
    new_sed = channel.apply_spectral_response(ais.source_sed, obj_wavelength)
    new_sky_sed = channel.apply_spectral_response(ais.sky_sed, obj_wavelength)
    ais.apply_sparc4_spectral_response('photometric')
    assert np.allclose(ais.source_sed, new_sed)
    assert np.allclose(ais.sky_sed, new_sky_sed)


def test_apply_sparc4_spectral_response_polarimetric(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    ais.create_sky_sed(moon_phase)
    channel = Channel(channel_id)
    channel.write_sparc4_operation_mode('polarimetric', 'polarizer', 'quarter')
    new_sed = channel.apply_spectral_response(ais.source_sed, obj_wavelength)
    new_sky_sed = channel.apply_spectral_response(ais.sky_sed, obj_wavelength)
    ais.apply_sparc4_spectral_response('polarimetric', 'polarizer', 'quarter')
    assert np.allclose(ais.source_sed, new_sed)
    assert np.allclose(ais.sky_sed, new_sky_sed)


# # --------------------------------------------------------------------------


def test_integrate_sed(ais):
    ais.create_source_sed(calculation_method, magnitude,
                          wavelegnth_interval, temperature)
    ais.create_sky_sed(moon_phase)
    src = Source()
    wv, sed2 = src.calculate_sed(
        calculation_method, magnitude, wavelegnth_interval, temperature)
    star_photons_per_second = np.trapz(sed2, wv)
    ais._integrate_sed()
    sky = Sky()
    sky_sed2 = sky.calculate_sed(moon_phase, obj_wavelength)
    sky_photons_per_second = np.trapz(sky_sed2, obj_wavelength)
    assert np.allclose(ais.star_photons_per_second, star_photons_per_second)
    assert np.allclose(ais.sky_photons_per_second, sky_photons_per_second)


image_path = os.path.join('tests', 'fits')


def test_find_image_index(ais):
    index = ais._find_image_index(image_path)
    assert index == 3


def test_create_image_name(ais):
    image_name = ais._create_image_name(image_path)
    now = datetime.datetime.now()
    datetime_str = now.strftime('%Y%m%d_s4c1_000004.fits')
    assert image_name == datetime_str

# --------------------------- test create artificial image ----------------


# def test_create_artificial_image_phot():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.apply_sparc4_spectral_response(sparc4_operation_mode)
#     ais.create_artificial_image()


# def test_create_artificial_image_pol():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     sparc4_operation_mode = {
#         "acquisition_mode": "polarimetric",
#         "calibration_wheel": "empty",
#         "retarder": "quarter",
#     }
#     ais = Artificial_Image_Simulator(
#         ccd_operation_mode,
#         image_dir=os.path.join("FITS"),
#     )
#     ais.apply_sparc4_spectral_response(sparc4_operation_mode)
#     ais.create_artificial_image()


# def test_create_background_image():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.create_background_image()


# def test_creat_bias_image():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.create_bias_image()


# def test_creat_dark_image():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.create_dark_image()


# def test_creat_random_image():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.create_random_image(n=2)


# def test_creat_flat_image():
#     ccd_operation_mode = {
#         "em_mode": "Conv",
#         "em_gain": 1,
#         "preamp": 1,
#         "hss": 1,
#         "binn": 1,
#         "t_exp": 1,
#         "ccd_temp": -70,
#         "image_size": 100,
#     }
#     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
#     ais.create_flat_image()
