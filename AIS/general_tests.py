# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""

# import os
# from random import gauss

# -------------------------------------------------------------------------------------------------------
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd

# import Header
# from Artificial_Image_Simulator import Artificial_Image_Simulator
# from Background_Image import Background_Image
# from Channel_Creator import Concrete_Channel_1
# from SPARC4_Spectral_Response import Concrete_SPARC4_Spectral_Response_1
# from Spectrum_Calculation import Spectrum_Calculation
# from Telescope_Spectral_Response import Telescope_Spectral_Response
# from tests.SPARC4_SR_curves import *
# from tests.test_AIS_operation import multiply_matrices

# from tests.SPARC4_SR_curves import *
# from tests.test_AIS_operation import multiply_matrices
# from tests.test_S4_SR_class import c1_s4_sr

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


# l_init, l_final, l_step = 400, 1150, 50
# magnitude = 22
# wavelength_interval = range(l_init, l_final, l_step)
# n = len(wavelength_interval)
# sc = Spectrum_Calculation(5700, l_init, l_final, l_step)
# specific_flux = sc.calculate_specific_flux(magnitude)

# c1_s4sr = Concrete_SPARC4_Spectral_Response_1(wavelength_interval)
# c1_s4sr.write_specific_flux(specific_flux.copy())


# c1_s4sr.apply_calibration_wheel()
# specific_flux = multiply_matrices(calibration_wheel, specific_flux)

# c1_s4sr.apply_retarder()
# specific_flux = multiply_matrices(retarder, specific_flux)

# c1_s4sr.apply_analyser()
# specific_flux = multiply_matrices(analyser_ordinary_ray, specific_flux)

# c1_s4sr.apply_collimator()
# specific_flux = np.multiply(colimator_transmitance, specific_flux)

# c1_s4sr.apply_dichroic()
# specific_flux = np.multiply(dichroic_c1, specific_flux)

# c1_s4sr.apply_camera()
# specific_flux - np.multiply(camera_c1, specific_flux)

# c1_s4sr.apply_ccd()
# specific_flux = np.multiply(ccd_transmitance_c1, specific_flux)
# print(specific_flux[0, :], "\n", c1_s4sr.specific_ordinary_ray[0, :])
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


# l_init, l_final, l_step = 400, 1150, 50
# num = int((l_final - l_init) / l_step)
# wavelength_interval = np.linspace(l_init, l_final, num)
# sc = Spectrum_Calculation(5700, l_init, l_final, l_step)
# pps = sc.calculate_specific_flux(22)
# plt.plot(wavelength_interval, pps[0, :])
# plt.show()

# ------------------------------------------------------------------------------------------------------

# import astropy.io.fits as fits

# from Header import Header

# dic = {
#     "em_mode": 0,
#     "em_gain": 1,
#     "preamp": 1,
#     "hss": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "ccd_temp": -70,
#     "image_size": 100,
# }

# hdr = Header(dic, 3, 9914)
# hdr._read_spreadsheet()
# header = hdr.create_header()
# fits.writeto("image.fits", np.zeros((10, 10)), header=header, overwrite=True)


# -------------------------------------------------------------------------------

# import os

# from Artificial_Image_Simulator import Artificial_Image_Simulator

# dic = {
#     "em_mode": 0,
#     "em_gain": 1,
#     "preamp": 2,
#     "hss": 0.1,
#     "binn": 1,
#     "t_exp": 1,
#     "ccd_temp": -70,
#     "image_size": 1024,
# }


# ais = Artificial_Image_Simulator(
#     ccd_operation_mode=dic,
#     channel=1,
#     star_coordinates=[512, 512],
#     bias_level=500,
#     sparc4_operation_mode="phot",
#     image_dir=os.path.join("..", "FITS"),
#     wavelength_interval=(400, 1100, 50),
#     star_temperature=5700,
# )

# ais.create_artificial_image()


import os

import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import splev, splrep

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Atmosphere_Spectral_Response import Atmosphere_Spectral_Response
from AIS.Channel_Creator import Concrete_Channel_1
from AIS.SPARC4_Spectral_Response import (
    Concrete_SPARC4_Spectral_Response_1,
    sparc4_spectral_response,
)
from tests.AIS_spectral_response_curves import star_specific_flux, wavelength_interval

# -------------------------------------------------------------------------------

s4_sr = Concrete_Channel_1("pol", wavelength_interval)
new_flux = s4_sr.apply_sparc4_spectral_response(star_specific_flux)


# -------------------------------------------------------------------------------
# from Read_Noise_Calculation import Read_Noise_Calculation

# dic = {
#     "em_mode": 1,
#     "em_gain": 2,
#     "preamp": 1,
#     "hss": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "ccd_temp": -70,
#     "image_size": 1024,
# }

# rn_calc = Read_Noise_Calculation(dic, channel=4)
# rn_calc.calculate_read_noise()

# -------------------------------------------------------------------------------


# air_mass = 1

# spreadsheet_path = os.path.join(
#     "Atmosphere_Spectral_Response", "atmosphere_spectral_response.csv"
# )
# spreadsheet = pd.read_csv(spreadsheet_path)
# atm_wavelength_interval = [float(value) for value in spreadsheet["Wavelength"][1:]]
# photometric_extinction_coef = [
#     float(value) / 100 for value in spreadsheet["photometric"][1:]
# ]
# transmitance = [10 ** (-0.4 * k * air_mass) for k in photometric_extinction_coef]
# spl = splrep(atm_wavelength_interval, transmitance)
# atm_transmitance = splev(wavelength_interval, spl)
# plt.plot(wavelength_interval, atm_transmitance)
# plt.show()