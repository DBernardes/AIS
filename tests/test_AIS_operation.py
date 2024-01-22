# # -*- coding: utf-8 -*-
# """AIS operation tests.

# This script presents the tests of the AIS operation. These tests are:

#     - Create Image Name: tess of the function that creates the image name based
#     on the provided operation mode of the CCD

#     - Configura Gain: tess of the function that sets the CCD gain based
#     on the provided operation mode of the CCD


# Created on Fri Apr 16 09:10:51 2021

# @author: denis
# """


# import os
# import datetime
# import numpy as np
# import pytest
# from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
# from AIS.Spectral_Energy_Distribution import Source, Sky
# from AIS.Spectral_Response import Atmosphere, Telescope, Channel
# from copy import copy
# from AIS.Noise import Noise
# from AIS.Point_Spread_Function import Point_Spread_Function

# ccd_operation_mode = {
#     "em_mode": "Conv",
#     "em_gain": 1,
#     "preamp": 1,
#     "readout": 1,
#     "binn": 1,
#     "t_exp": 1,
#     "image_size": 100,
# }

# calculation_method = "blackbody"
# magnitude = 15
# wavelegnth_interval = (400, 1100, 100)
# star_temperature = 5700
# ccd_temperature = -70
# channel_id = 1


# @pytest.fixture
# def ais():
#     return Artificial_Image_Simulator(ccd_operation_mode, channel_id, ccd_temperature)


# # ------------------------------------------------------------


# def test_print_available_spectral_types(ais):
#     ais.print_available_spectral_types()


# def test_create_source_sed_blackbody(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     src = Source()
#     wv2, sed2 = src.calculate_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     assert np.allclose(ais.wavelength, wv2)
#     assert np.allclose(ais.source_sed, sed2)


# def test_create_source_sed_spectral_lib(ais):
#     calculation_method = "spectral_library"
#     spectral_type = "A0V"
#     ais.create_source_sed(calculation_method, magnitude, spectral_type=spectral_type)
#     src = Source()
#     wv2, sed2 = src.calculate_sed(
#         calculation_method, magnitude, spectral_type=spectral_type
#     )
#     assert np.allclose(ais.wavelength, wv2)
#     assert np.allclose(ais.source_sed, sed2)


# def test_write_source_sed(ais):
#     n = len(wavelegnth_interval)
#     sed = np.ones((4, n))
#     ais.write_source_sed(wavelegnth_interval, sed)
#     assert np.allclose(wavelegnth_interval, ais.wavelength)
#     assert np.allclose(sed, ais.source_sed)


# wv = wavelegnth_interval
# obj_wavelength = np.linspace(wv[0], wv[1], wv[2])
# moon_phase = "new"


# def test_create_sky_sed(ais):
#     ais.create_source_sed("blackbody", magnitude, wavelegnth_interval, star_temperature)
#     ais.create_sky_sed(moon_phase)
#     sky = Sky()
#     sed2 = sky.calculate_sed(moon_phase, obj_wavelength)

#     assert np.allclose(ais.sky_sed, sed2)


# # # ----------------------- Apply spectruns ---------------------------
# air_mass = 1
# sky_condition = "photometric"


# def test_apply_atmosphere_spectral_response(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     atm = Atmosphere()
#     new_sed = atm.apply_spectral_response(
#         obj_wavelength, copy(ais.source_sed), air_mass, sky_condition
#     )

#     ais.apply_atmosphere_spectral_response(air_mass, sky_condition)
#     assert np.allclose(ais.source_sed, new_sed)


# def test_apply_telescope_spectral_response(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)
#     tel = Telescope()
#     new_sed = tel.apply_spectral_response(obj_wavelength, ais.source_sed)
#     new_sky_sed = tel.apply_spectral_response(obj_wavelength, ais.sky_sed)
#     ais.apply_telescope_spectral_response()
#     assert np.allclose(ais.source_sed, new_sed)
#     assert np.allclose(ais.sky_sed, new_sky_sed)


# # -----------------------------------Test apply SPARC4 spectral response ----------------------------------------------------


# def test_apply_sparc4_spectral_response_photometric(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)
#     channel = Channel(channel_id)
#     channel.write_sparc4_operation_mode("photometry")
#     new_sed = channel.apply_spectral_response(ais.source_sed, obj_wavelength)
#     new_sky_sed = channel.apply_spectral_response(ais.sky_sed, obj_wavelength)
#     ais.apply_sparc4_spectral_response("photometry")
#     assert np.allclose(ais.source_sed, new_sed)
#     assert np.allclose(ais.sky_sed, new_sky_sed)


# def test_apply_sparc4_spectral_response_polarimetric(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)

#     channel = Channel(channel_id)
#     channel.write_sparc4_operation_mode("polarimetry", "polarizer", "quarter")
#     new_sed = channel.apply_spectral_response(copy(ais.source_sed), obj_wavelength)
#     new_sky_sed = channel.apply_spectral_response(copy(ais.sky_sed), obj_wavelength)

#     ais.apply_sparc4_spectral_response("polarimetry", "polarizer", "quarter")

#     assert np.allclose(ais.source_sed, new_sed)
#     assert np.allclose(ais.sky_sed, new_sky_sed)


# # # --------------------------------------------------------------------------


# def test_integrate_sed(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)
#     src = Source()
#     wv, sed2 = src.calculate_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     star_photons_per_second = np.trapz(sed2, wv * 1e-9)
#     ais._integrate_sed()
#     sky = Sky()
#     sky_sed2 = sky.calculate_sed(moon_phase, obj_wavelength)
#     sky_photons_per_second = np.trapz(sky_sed2, obj_wavelength * 1e-9)
#     assert np.allclose(ais.star_photons_per_second, star_photons_per_second)
#     assert np.allclose(ais.sky_photons_per_second, sky_photons_per_second)


# image_path = os.path.join("tests", "fits")


# def test_find_image_index(ais):
#     index = ais._find_image_index(image_path)
#     assert index == 3


# def test_create_image_name(ais):
#     ais._create_image_name(image_path)
#     now = datetime.datetime.now()
#     datetime_str = now.strftime("%Y%m%d_s4c1_000004.fits")
#     assert ais.image_name == datetime_str


# # --------------------------- test create artificial image ----------------


# sparc4_operation_mode = "photometry"
# star_coordinates = (50, 50)

# path = os.path.join("tests", "fits")


# def test_create_artificial_image_phot(ais):
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)
#     ais.create_artificial_image(path, star_coordinates)
#     os.remove(os.path.join(path, ais.image_name))


# # def test_create_artificial_image_pol():
# #     ccd_operation_mode = {
# #         "em_mode": "Conv",
# #         "em_gain": 1,
# #         "preamp": 1,
# #         "hss": 1,
# #         "binn": 1,
# #         "t_exp": 1,
# #         "ccd_temp": -70,
# #         "image_size": 100,
# #     }
# #     sparc4_operation_mode = {
# #         "acquisition_mode": "polarimetric",
# #         "calibration_wheel": "empty",
# #         "retarder": "quarter",
# #     }
# #     ais = Artificial_Image_Simulator(
# #         ccd_operation_mode,
# #         image_dir=os.path.join("FITS"),
# #     )
# #     ais.apply_sparc4_spectral_response(sparc4_operation_mode)
# #     ais.create_artificial_image()


# def test_create_background_image():
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     ais.create_source_sed(
#         calculation_method, magnitude, wavelegnth_interval, star_temperature
#     )
#     ais.create_sky_sed(moon_phase)
#     ais.apply_sparc4_spectral_response("photometry")
#     ais.create_background_image(path)
#     os.remove(os.path.join(path, ais.image_name))


# def test_creat_background_image_error(ais):
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     with pytest.raises(ValueError):
#         ais.create_background_image(path, 0)


# def test_creat_bias_image(ais):
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     ais.create_bias_image(path)
#     os.remove(os.path.join(path, ais.image_name))


# def test_creat_bias_image_error(ais):
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     with pytest.raises(ValueError):
#         ais.create_bias_image(path, 0)


# def test_creat_dark_image():
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     ais.create_dark_image(path)
#     os.remove(os.path.join(path, ais.image_name))


# def test_creat_dark_image_error(ais):
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     with pytest.raises(ValueError):
#         ais.create_dark_image(path, 0)


# # def test_creat_random_image():
# #     ccd_operation_mode = {
# #         "em_mode": "Conv",
# #         "em_gain": 1,
# #         "preamp": 1,
# #         "hss": 1,
# #         "binn": 1,
# #         "t_exp": 1,
# #         "ccd_temp": -70,
# #         "image_size": 100,
# #     }
# #     ais = Artificial_Image_Simulator(ccd_operation_mode, image_dir=os.path.join("FITS"))
# #     ais.create_random_image(n=2)


# def test_creat_flat_image():
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     ais.create_flat_image(path)
#     os.remove(os.path.join(path, ais.image_name))


# def test_creat_flat_image_error(ais):
#     ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
#     with pytest.raises(ValueError):
#         ais.create_flat_image(path, 0)


# def test_solve_second_degree_equation(ais):
#     roots = ais._solve_second_degree_eq(1, -4, 4)
#     assert roots == [2, 2]


# def test_calculate_exposure_time():
#     ais = Artificial_Image_Simulator(
#         ccd_operation_mode, channel_id=1, ccd_temperature=-70
#     )
#     ais.create_source_sed(
#         calculation_method="spectral_library",
#         magnitude=15,
#         wavelength_interval=(400, 1100, 100),
#         spectral_type="A0v",
#     )
#     ais.create_sky_sed(moon_phase="new")
#     ais.apply_atmosphere_spectral_response()
#     ais.apply_telescope_spectral_response()
#     ais.apply_sparc4_spectral_response(acquisition_mode="photometry")
#     ais_texp = ais.calculate_exposure_time()

#     snr = 100
#     psf_obj = Point_Spread_Function(ccd_operation_mode, channel_id)
#     n_pix = psf_obj.calculate_npix_star(seeing=1.5)
#     noise_obj = Noise(channel_id)
#     dark_noise = noise_obj.calculate_dark_current(-70)
#     read_noise = noise_obj.calculate_read_noise(ccd_operation_mode)
#     noise_factor = 1
#     em_gain = 1
#     binn = ccd_operation_mode["binn"]
#     if ccd_operation_mode["em_mode"] == "EM":
#         noise_factor = 1.4
#         em_gain = ccd_operation_mode["em_gain"]
#     star_photons_per_second = ais.star_photons_per_second
#     sky_photons_per_second = ais.sky_photons_per_second

#     a = star_photons_per_second**2
#     b = (
#         snr**2
#         * noise_factor**2
#         * (star_photons_per_second + n_pix * (sky_photons_per_second + dark_noise))
#     )
#     c = snr**2 * n_pix * (read_noise / em_gain / binn) ** 2

#     texp_list = ais._solve_second_degree_eq(a, -b, -c)
#     min_t_exp = min([texp for texp in texp_list if texp > 0])

#     assert min_t_exp == ais_texp
