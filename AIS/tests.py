# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

import os
from random import gauss

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import Header
from Artificial_Image_Simulator import Artificial_Image_Simulator
from Background_Image import Background_Image
from Channel_Creator import Concrete_Channel_1
from SPARC4_Spectral_Response import Concrete_SPARC4_Spectral_Response_1
from Spectrum_Calculation import Spectrum_Calculation
from Telescope_Spectral_Response import Telescope_Spectral_Response
from tests.SPARC4_SR_curves import *
from tests.test_AIS_operation import multiply_matrices
from tests.test_S4_SR_class import c1_s4_sr

# dic = {
#     "em_mode": 0,
#     "em_gain": 1,
#     "preamp": 1,
#     "hss": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "ccd_temp": -70,
#     "image_size": 200,
# }

# gaussian_std = 8
# ais = Artificial_Image_Simulator(
#     dic,
#     image_dir=os.path.join("..", "FITS"),
#     sparc4_operation_mode="pol",
#     gaussian_std=gaussian_std,
# )
# ais.apply_sparc4_spectral_response()

# -------------------------------------------------------------------------------

l_init, l_final, l_step = 400, 1150, 50
wavelength_interval = range(l_init, l_final, l_step)
n = len(wavelength_interval)
sc = Spectrum_Calculation(5700, l_init, l_final, l_step)
specific_flux = sc.calculate_star_specific_flux()

# c1_s4sr = Concrete_SPARC4_Spectral_Response_1()
# c1_s4sr.write_specific_flux(specific_flux.copy(), wavelength_interval)

# c1_s4sr.apply_calibration_wheel()
# specific_flux = multiply_matrices(calibration_wheel, specific_flux)

# c1_s4sr.apply_retarder()
# specific_flux = multiply_matrices(retarder, specific_flux)

# c1_s4sr.apply_analyser()
# specific_flux = multiply_matrices(analyser_ordinary_ray, specific_flux)

# c1_s4sr.apply_collimator()
# specific_flux = np.multiply(colimator_transmitance, specific_flux[0, :])

# c1_s4sr.apply_dichroic()
# specific_flux = np.multiply(dichroic_c1, specific_flux)

# c1_s4sr.apply_camera()
# specific_flux - np.multiply(camera_c1, specific_flux)

# c1_s4sr.apply_ccd()
# specific_flux = np.multiply(ccd_transmitance_c1, specific_flux)
# print(np.allclose(specific_flux, c1_s4sr.specific_ordinary_ray))

# -------------------------------------------------------------------------------


# c1_extra_ordinary_ray_phot = 0


# c1_s4sr = Concrete_SPARC4_Spectral_Response_1()
# c1_s4sr.write_specific_flux(specific_flux.copy(), wavelength_interval)

# c1_s4sr.apply_collimator()
# c1_ordinary_ray_phot = np.multiply(colimator_transmitance, specific_flux.copy()[0, :])


# c1_s4sr.apply_dichroic()
# c1_ordinary_ray_phot = np.multiply(dichroic_c1, c1_ordinary_ray_phot)


# c1_s4sr.apply_camera()
# c1_ordinary_ray_phot = np.multiply(camera_c1, c1_ordinary_ray_phot)


# c1_s4sr.apply_ccd()
# c1_ordinary_ray_phot = np.multiply(ccd_transmitance_c1, c1_ordinary_ray_phot)
# print(np.allclose(c1_s4sr.specific_ordinary_ray, c1_ordinary_ray_phot))


# chc1 = Concrete_Channel_1("phot")
# (
#     specific_star_ordinary_ray,
#     specific_star_extra_ordinary_ray,
# ) = chc1.apply_sparc4_spectral_response(np.ones((4, 16)), 350, 1150, 50)
# print(specific_star_ordinary_ray)

# -------------------------------------------------------------------------------

# ccd_operation_mode = {
#     "em_mode": 0,
#     "em_gain": 1,
#     "preamp": 1,
#     "hss": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "image_size": 1024,
# }

# sky_flux = 10
# ccd_gain = 3
# dark_current = 1e-5
# read_noise = 6.67
# bias_level = 500
# dc = dark_current * 1
# rn = read_noise
# nf = 1

# noise = (
#     np.sqrt(
#         rn ** 2
#         + (sky_flux + dc)
#         * ccd_operation_mode["t_exp"]
#         * nf ** 2
#         * ccd_operation_mode["em_gain"] ** 2
#         * ccd_operation_mode["binn"] ** 2
#     )
#     / ccd_gain
# )

# bgi = Background_Image(
#     ccd_operation_mode, ccd_gain, dark_current, read_noise, bias_level
# )
# image = bgi.create_background_image(sky_flux)
# mean = np.mean(image)
# std = np.std(image)
# print(std, noise)
# print(np.allclose(std, noise, atol=0.005))

# -------------------------------------------------------------------------------------------------------

# sc = Spectrum_Calculation(5700, l_init, l_final, l_step)
# pps = sc.calculate_star_specific_flux()
# print(pps[0, :])

# ------------------------------------------------------------------------------------------------------

import astropy.io.fits as fits

from Header import Header

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

hdr = Header(dic, 3, 9914)
hdr._read_spreadsheet()
header = hdr.create_header()
fits.writeto("image.fits", np.zeros((10, 10)), header=header, overwrite=True)
