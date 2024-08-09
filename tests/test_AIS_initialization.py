# # -*- coding: utf-8 -*-
"""AIS operation tests.

This script presents the tests of the AIS initialization. These tests are:

    - Initialization of the AIS class: it is tested if the the values of the
    provided parameters are correct. The allowed parameter are the star flux,
    the sky flux, the gaussian distribution of the PSF, the bias level of the
    image, and the image directory where the image should be saved. The value
    of the numerical parameters should be greater than zero. The value of the
    image directory should be a string.

    - Operation Mode of the CCD: the parameters of the operation mode of the
    CCD are provided as a python ccd_operation_modetionary. It is tested if the keyword values
    of the provided ccd_operation_modetionary are correct; it is tested if the values of the
    parameters are in agreement with the manufacture manual.

    - Channel ID: each SPARC4 channel is represented as a python class, and it
    receives a channel identifier (ID). This test calls the respective object
    for the channel provided by the user (1, 2, 3, or 4) and executs its
    function get_channel_ID(self). The returned results should be in agreement with
    the provided channel.


Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


import datetime
import os
import unittest

import astropy.io.fits as fits
import numpy as np
import pytest

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Background_Image import Background_Image
from AIS.Header import Header
from AIS.Noise import Noise
from AIS.Point_Spread_Function import Point_Spread_Function


class Test_AIS_Initialization(unittest.TestCase):
    CH_ID = 1
    TEMP = -70
    CCD_OP_MODE = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "readout": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 200,
    }

    AIS = Artificial_Image_Simulator(CCD_OP_MODE, CH_ID, TEMP)

    @classmethod
    def setUpClass(cls):
        pass

    # -------------------- Testing the AIS class --------------------------------

    def test_ccd_operation_mode(self):
        assert self.AIS.ccd_operation_mode == self.CCD_OP_MODE

    def test_channel_id(self):
        assert self.AIS.channel_id == self.CH_ID

    # -------  provide a wrong parameter to the CCD operation mode---------------

    def test_ccd_operation_mode_wrong_keyword_em_mode(self):
        self.CCD_OP_MODE["em_modee"] = "Conv"
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_em_gain(self):
        self.CCD_OP_MODE["em_gainn"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_preamp(self):
        self.CCD_OP_MODE["preampp"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_hss(self):
        self.CCD_OP_MODE["readout"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_bin(self):
        self.CCD_OP_MODE["bin"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_t_exp(self):
        self.CCD_OP_MODE["texpp"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_ccd_temp(self):
        self.CCD_OP_MODE["ccd_tempp"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyword_image_size(self):
        self.CCD_OP_MODE["image_sizee"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    # ------------ provide a wrong value to the CCD operation mode-------------

    def test_ccd_operation_mode_wrong_keyvalue_em_mode(self):
        self.CCD_OP_MODE["em_mode"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_em_gain_1(self):
        self.CCD_OP_MODE["em_gain"] = 0
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_em_gain_2(self):
        self.CCD_OP_MODE["em_gain"] = 500
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_em_gain_3(self):
        self.CCD_OP_MODE["em_gain"] = 2
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_preamp(self):
        self.CCD_OP_MODE["preamp"] = 4
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_hss(self):
        self.CCD_OP_MODE["hss"] = -1
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_bin(self):
        self.CCD_OP_MODE["binn"] = 0
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_t_exp_1(self):
        self.CCD_OP_MODE["t_exp"] = 0
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_wrong_keyvalue_ccd_temp_1(self):
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, 0)

    def test_ccd_operation_mode_wrong_keyvalue_t_image_size_1(self):
        self.CCD_OP_MODE["image_size"] = 0
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    # ---------------------------Missing CCD operation mode parameter -------------

    def test_ccd_operation_mode_missing_parameter_em_mode(self):
        del self.CCD_OP_MODE["em_mode"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_em_gain(self):
        del self.CCD_OP_MODE["em_gain"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_preamp(self):
        del self.CCD_OP_MODE["preamp"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_hss(self):
        del self.CCD_OP_MODE["readout"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_t_exp(self):
        del self.CCD_OP_MODE["t_exp"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_bin(self):
        del self.CCD_OP_MODE["binn"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)

    def test_ccd_operation_mode_missing_parameter_image_size(self):
        del self.CCD_OP_MODE["image_size"]
        with pytest.raises(ValueError):
            Artificial_Image_Simulator(self.CCD_OP_MODE, self.CH_ID, self.TEMP)


class Test_AIS_Operation(unittest.TestCase):
    CH_ID = 1
    TEMP = -70
    CCD_OP_MODE = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "readout": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 200,
    }
    AIS = Artificial_Image_Simulator(CCD_OP_MODE, CH_ID, TEMP)
    MAGNITUDE = 15
    FITS_PATH = os.path.join("tests", "fits")
    SEED = 5

    @classmethod
    def setUpClass(cls):
        pass

    def test_print_available_spectral_types(self):
        self.AIS.print_available_spectral_types()

    def test_verify_var_interval(self):
        self.AIS._verify_var_in_interval(100, "dummy", 0, 200)

    def test_verify_var_interval_2(self):
        with pytest.raises(ValueError):
            self.AIS._verify_var_in_interval(101, "dummy", 0, 100)

    def test_verify_var_interval_3(self):
        with pytest.raises(ValueError):
            self.AIS._verify_var_in_interval(-1, "dummy", 0, 100)

    def test_check_var_in_a_list(self):
        self.AIS._check_var_in_a_list("a", "dummy", ["a"])

    def test_check_var_in_a_list_2(self):
        with pytest.raises(ValueError):
            self.AIS._check_var_in_a_list("a", "dummy", [])

    def test_write_source_sed(self):
        wavelength = np.linspace(400, 1100, 100)
        sed = np.ones((4, 100))
        self.AIS.write_source_sed(wavelength, sed)

        assert np.allclose(sed, self.AIS.source_sed)

    def test_write_source_sed_1(self):
        wavelength = np.linspace(400, 1100, 100)
        sed = np.zeros((100))
        self.AIS.write_source_sed(wavelength, sed)

        assert np.allclose(np.zeros((4, 100)), self.AIS.source_sed)

    def test_integrate_sed(self):
        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")

        star_photons_per_second = np.trapezoid(
            self.AIS.source_sed, self.AIS.wavelength * 1e-9
        )
        sky_photons_per_second = np.trapezoid(
            self.AIS.sky_sed, self.AIS.wavelength * 1e-9
        )
        self.AIS._integrate_sed()

        assert np.allclose(self.AIS.star_photons_per_second, star_photons_per_second)
        assert np.allclose(self.AIS.sky_photons_per_second, sky_photons_per_second)

    def test_find_image_index(self):
        index = self.AIS._find_image_index(self.FITS_PATH)
        assert index == 3

    def test_create_image_name(self):
        self.AIS._create_image_name(self.FITS_PATH)
        now = datetime.datetime.now()
        datetime_str = now.strftime("%Y%m%d_s4c1_000004.fits")
        assert self.AIS.image_name == datetime_str

    def test_calculate_exposure_time(self):
        snr = 100
        seeing = 1
        psf = Point_Spread_Function(self.CCD_OP_MODE, self.CH_ID)
        n_pix = psf.calculate_npix_star(seeing)
        noise_obj = Noise(self.CH_ID)
        dark_noise = noise_obj.calculate_dark_current(self.TEMP)
        read_noise = noise_obj.calculate_read_noise(self.CCD_OP_MODE)
        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("photometry")
        min_t_exp = self.AIS.calculate_exposure_time(snr, seeing)

        noise_factor = 1.0
        em_gain = 1.0
        binn = self.CCD_OP_MODE["binn"]

        a = self.AIS.star_photons_per_second**2
        b = (
            snr**2
            * noise_factor**2
            * (
                self.AIS.star_photons_per_second
                + n_pix * (self.AIS.sky_photons_per_second + dark_noise)
            )
        )
        c = snr**2 * n_pix * (read_noise / em_gain / binn) ** 2

        texp_list = self.AIS._solve_second_degree_eq(a, -b, -c)

        assert min_t_exp == min([texp for texp in texp_list if texp > 0])

    # def test_create_source_sed_blackbody(self):
    #     self.AIS.create_source_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     src = Source()
    #     wv2, sed2 = src.calculate_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     assert np.allclose(self.AIS.wavelength, wv2)
    #     assert np.allclose(self.AIS.source_sed, sed2)

    # def test_create_source_sed_spectral_lib(self):
    #     calculation_method = "spectral_library"
    #     spectral_type = "A0V"
    #     self.AIS.create_source_sed(calculation_method, magnitude, spectral_type=spectral_type)
    #     src = Source()
    #     wv2, sed2 = src.calculate_sed(
    #         calculation_method, magnitude, spectral_type=spectral_type
    #     )
    #     assert np.allclose(self.AIS.wavelength, wv2)
    #     assert np.allclose(self.AIS.source_sed, sed2)

    # def test_create_sky_sed(self):
    #     self.AIS.create_source_sed("blackbody", magnitude, wavelegnth_interval, star_temperature)
    #     self.AIS.create_sky_sed(moon_phase)
    #     sky = Sky()
    #     sed2 = sky.calculate_sed(moon_phase, obj_wavelength)

    #     assert np.allclose(self.AIS.sky_sed, sed2)

    # # # ----------------------- Apply spectruns ---------------------------
    # air_mass = 1
    # sky_condition = "photometric"

    # def test_apply_atmosphere_spectral_response(self):
    #     self.AIS.create_source_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     atm = Atmosphere()
    #     new_sed = atm.apply_spectral_response(
    #         obj_wavelength, copy(self.AIS.source_sed), air_mass, sky_condition
    #     )

    #     self.AIS.apply_atmosphere_spectral_response(air_mass, sky_condition)
    #     assert np.allclose(self.AIS.source_sed, new_sed)

    # def test_apply_telescope_spectral_response(self):
    #     self.AIS.create_source_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     self.AIS.create_sky_sed(moon_phase)
    #     tel = Telescope()
    #     new_sed = tel.apply_spectral_response(obj_wavelength, self.AIS.source_sed)
    #     new_sky_sed = tel.apply_spectral_response(obj_wavelength, self.AIS.sky_sed)
    #     self.AIS.apply_telescope_spectral_response()
    #     assert np.allclose(self.AIS.source_sed, new_sed)
    #     assert np.allclose(self.AIS.sky_sed, new_sky_sed)

    # # -----------------------------------Test apply SPARC4 spectral response ----------------------------------------------------

    # def test_apply_sparc4_spectral_response_photometric(self):
    #     self.AIS.create_source_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     self.AIS.create_sky_sed(moon_phase)
    #     channel = Channel(channel_id)
    #     channel.write_sparc4_operation_mode("photometry")
    #     new_sed = channel.apply_spectral_response(self.AIS.source_sed, obj_wavelength)
    #     new_sky_sed = channel.apply_spectral_response(self.AIS.sky_sed, obj_wavelength)
    #     self.AIS.apply_sparc4_spectral_response("photometry")
    #     assert np.allclose(self.AIS.source_sed, new_sed)
    #     assert np.allclose(self.AIS.sky_sed, new_sky_sed)

    # def test_apply_sparc4_spectral_response_polarimetric(self):
    #     self.AIS.create_source_sed(
    #         calculation_method, magnitude, wavelegnth_interval, star_temperature
    #     )
    #     self.AIS.create_sky_sed(moon_phase)

    #     channel = Channel(channel_id)
    #     channel.write_sparc4_operation_mode("polarimetry", "polarizer", "quarter")
    #     new_sed = channel.apply_spectral_response(copy(self.AIS.source_sed), obj_wavelength)
    #     new_sky_sed = channel.apply_spectral_response(copy(self.AIS.sky_sed), obj_wavelength)

    #     self.AIS.apply_sparc4_spectral_response("polarimetry", "polarizer", "quarter")

    #     assert np.allclose(self.AIS.source_sed, new_sed)
    #     assert np.allclose(self.AIS.sky_sed, new_sky_sed)

    # # --------------------------- test create artificial image ----------------

    def test_create_artificial_image_phot(self):
        star_coordinates = (50, 50)
        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("photometry")
        self.AIS.create_artificial_image(self.FITS_PATH, star_coordinates, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("photometry")
        self.AIS._integrate_sed()
        ord_ray, extra_ord_ray = self.AIS.star_photons_per_second, 0

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        background = bgi.create_sky_background(
            self.AIS.sky_photons_per_second, self.SEED
        )
        psf = Point_Spread_Function(self.CCD_OP_MODE, self.CH_ID)
        star_psf = psf.create_star_image(
            star_coordinates, ord_ray, extra_ord_ray, 1, self.SEED
        )
        self.AIS._create_image_name(self.FITS_PATH)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        header["OBSTYPE"] = "OBJECT"
        header["FILENAME"] = self.AIS.image_name
        header["SHUTTER"] = "OPEN"

        image = background + star_psf
        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)
        fits.writeto(
            file2,
            image,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    def test_create_artificial_image_pol(self):
        star_coordinates = (50, 50)
        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("polarimetry")
        self.AIS.create_artificial_image(self.FITS_PATH, star_coordinates, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("polarimetry")
        self.AIS._integrate_sed()
        self.AIS.sky_photons_per_second = sum(self.AIS.sky_photons_per_second)
        ord_ray, extra_ord_ray = (
            self.AIS.star_photons_per_second[0],
            self.AIS.star_photons_per_second[1],
        )

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        background = bgi.create_sky_background(
            self.AIS.sky_photons_per_second, self.SEED
        )
        psf = Point_Spread_Function(self.CCD_OP_MODE, self.CH_ID)
        star_psf = psf.create_star_image(
            star_coordinates, ord_ray, extra_ord_ray, 1, self.SEED
        )
        self.AIS._create_image_name(self.FITS_PATH)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        header["OBSTYPE"] = "OBJECT"
        header["FILENAME"] = self.AIS.image_name
        header["SHUTTER"] = "OPEN"

        image = background + star_psf
        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)
        fits.writeto(
            file2,
            image,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    def test_create_background_image(self):
        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("photometry")
        self.AIS.create_background_image(self.FITS_PATH, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        self.AIS.create_source_sed("blackbody", self.MAGNITUDE, (400, 1100, 100), 5700)
        self.AIS.create_sky_sed("full")
        self.AIS.apply_sparc4_spectral_response("photometry")
        self.AIS._integrate_sed()

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        background = bgi.create_sky_background(
            self.AIS.sky_photons_per_second, self.SEED
        )

        self.AIS._create_image_name(self.FITS_PATH)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        header["OBSTYPE"] = "OBJECT"
        header["FILENAME"] = self.AIS.image_name
        header["SHUTTER"] = "OPEN"

        image = background
        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)
        fits.writeto(
            file2,
            image,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    def test_creat_bias_image(self):
        self.AIS.create_bias_image(self.FITS_PATH, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        bias = bgi.create_bias_background(self.SEED)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        self.AIS._create_image_name(self.FITS_PATH)
        header["EXPTIME"] = 1e-5
        header["FILENAME"] = self.AIS.image_name
        header["OBSTYPE"] = "ZERO"

        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        fits.writeto(
            file2,
            bias,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    def test_create_dark_image(self):
        self.AIS.create_dark_image(self.FITS_PATH, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        dark = bgi.create_dark_background(self.SEED)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        self.AIS._create_image_name(self.FITS_PATH)
        header["EXPTIME"] = 1e-5
        header["FILENAME"] = self.AIS.image_name
        header["OBSTYPE"] = "ZERO"

        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        fits.writeto(
            file2,
            dark,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    def test_creat_flat_image(self):
        self.AIS.create_flat_image(self.FITS_PATH, 1, self.SEED)
        file1 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        bgi = Background_Image(self.CCD_OP_MODE, self.CH_ID, self.TEMP)
        flat = bgi.create_flat_background(self.SEED)
        hdr = Header(self.CCD_OP_MODE, self.TEMP, self.CH_ID)
        header = hdr.create_header()
        self.AIS._create_image_name(self.FITS_PATH)
        header["EXPTIME"] = 1e-5
        header["FILENAME"] = self.AIS.image_name
        header["OBSTYPE"] = "ZERO"

        file2 = os.path.join(self.FITS_PATH, self.AIS.image_name)

        fits.writeto(
            file2,
            flat,
            header=header,
        )
        with fits.open(file1, memmap=True) as hdulist:
            image1 = hdulist[0].data.copy()
        with fits.open(file2, memmap=True) as hdulist:
            image2 = hdulist[0].data.copy()

        assert np.allclose(image1, image2)
        os.remove(file1)
        os.remove(file2)

    # def test_solve_second_degree_equation(self):
    #     roots = self.AIS._solve_second_degree_eq(1, -4, 4)
    #     assert roots == [2, 2]
