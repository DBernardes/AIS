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
from AIS.Spectral_Energy_Distribution import Source

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
    magnitude,
    multiply_matrices,
    polarizer_transmitance,
    retardance_quarter,
    retarder_transmitance,
    sky_condition,
    sky_specific_photons_per_second,
    sparc4_operation_mode,
    star_specific_photons_per_second,
    tel_reflectance,
)


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=ccd_operation_mode)


calculation_method = 'blackbody'
magnitude = 10
wavelegnth_interval = (350, 1100, 100)
temperature = 5700

# ------------------------------------------------------------


def test_create_source_sed_blackbody(ais):
    wv, sed = ais.create_source_sed(calculation_method, magnitude,
                                    wavelegnth_interval, temperature)
    src = Source()
    wv2, sed2 = src.calculate_sed(
        calculation_method, magnitude, wavelegnth_interval, temperature)
    assert np.allclose(wv, wv2)
    assert np.allclose(sed, sed2)


calculation_method = 'spectral_library'
spectral_type = 'O'


def test_create_source_sed_spectral_lib(ais):
    wv, sed = ais.create_source_sed(
        calculation_method, magnitude, spectral_type=spectral_type)
    src = Source()
    wv2, sed2 = src.calculate_sed(
        calculation_method, magnitude, spectral_type=spectral_type)
    assert np.allclose(wv, wv2)
    assert np.allclose(sed, sed2)


# # ----------------------- Apply spectruns ---------------------------


# def test_apply_atmosphere_spectral_response_star(ais):
#     new_star_specific_photons_per_second = np.multiply(
#         star_specific_photons_per_second[0][0].copy(), atm_transmitance
#     )
#     ais.apply_atmosphere_spectral_response(air_mass, sky_condition)
#     assert np.allclose(
#         ais.star_specific_photons_per_second[0][0],
#         new_star_specific_photons_per_second,
#     )


# def test_apply_atmosphere_spectral_response_sky(ais):
#     new_sky_specific_photons_per_second = np.multiply(
#         sky_specific_photons_per_second[0][0].copy(), atm_transmitance
#     )
#     ais.apply_atmosphere_spectral_response(air_mass, sky_condition)
#     assert np.allclose(
#         ais.sky_specific_photons_per_second[0][0].copy(),
#         new_sky_specific_photons_per_second,
#     )


# def test_apply_telescope_spectral_response_star(ais):
#     ais.apply_telescope_spectral_response()
#     new_star_specific_photons_per_second = np.multiply(
#         star_specific_photons_per_second.copy(), tel_reflectance
#     )
#     assert np.allclose(
#         ais.star_specific_photons_per_second, new_star_specific_photons_per_second
#     )


# def test_apply_telescope_spectral_response_sky(ais):
#     new_sky_specific_photons_per_second = np.multiply(
#         sky_specific_photons_per_second.copy(), tel_reflectance
#     )
#     ais.apply_telescope_spectral_response()
#     assert np.allclose(
#         ais.sky_specific_photons_per_second, new_sky_specific_photons_per_second
#     )


# # -----------------------------------Test apply SPARC4 spectral response ----------------------------------------------------


# # def test_apply_sparc4_spectral_response_star():
# #     star_specific_photons_per_second_c1 = star_specific_photons_per_second.copy()
# #     # Applying collimator
# #     star_specific_photons_per_second_c1[0][0] = np.multiply(
# #         collimator_transmitance, star_specific_photons_per_second_c1[0][0]
# #     )
# #     # Applying camera
# #     star_specific_photons_per_second_c1[0][0] = np.multiply(
# #         camera_c1, star_specific_photons_per_second_c1[0][0]
# #     )
# #     # Applying ccd
# #     star_specific_photons_per_second_c1[0][0] = np.multiply(
# #         ccd_transmitance_c1, star_specific_photons_per_second_c1[0][0]
# #     )

# #     sparc4_operation_mode = {
# #         "acquisition_mode": "photometric",
# #     }
# #     ais = Artificial_Image_Simulator(
# #         ccd_operation_mode=ccd_operation_mode,
# #         sparc4_operation_mode=sparc4_operation_mode,
# #         channel=1,
# #     )
# #     ais.apply_sparc4_spectral_response()
# #     assert np.allclose(
# #         star_specific_photons_per_second_c1, ais.star_specific_photons_per_second
# #     )


# # ----------------------------------------------------------------------------------------------------------------


# # def test_apply_sparc4_spectral_response_sky(ais):
# #     ais.apply_sparc4_spectral_response()

# #     assert np.allclose(ais.sky_specific_photons_per_second, new_ordinary_ray)


# # -----------------------------test _integrate_fluxes ------------------------


# def test_integrate_specific_photons_per_second(ais):
#     photons_per_second = np.trapz(star_specific_photons_per_second[0][0, :])
#     ais.star_specific_photons_per_second = star_specific_photons_per_second
#     ais.sky_specific_photons_per_second = star_specific_photons_per_second
#     ais._integrate_specific_photons_per_second()
#     assert ais.star_ordinary_ray == photons_per_second
#     assert ais.sky_photons_per_second == 25


# # --------------------------- test create artificial image ----------------


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
