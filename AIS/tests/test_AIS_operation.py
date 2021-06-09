# -*- coding: utf-8 -*-
"""AIS operation tests.

This script presents the tests of the AIS operation. These tests are:

    - Initialization of the AIS class: it is tested if the the values of the
    provided parameters are correct. The allowed parameter are the star flux,
    the sky flux, the gaussian distribution of the PSF, the bias level of the
    image, and the image directory where the image should be saved. The value
    of the numerical parameters should be greater than zero. The value of the
    image directory should be a string.

    - Operation Mode of the CCD: the parameters of the operation mode of the
    CCD are provided as a python dictionary. It is tested if the keyword values
    of the provided dictionary are correct; it is tested if the values of the
    parameters are in agreement with the manufacture manual.

    - Channel ID: each SPARC4 channel is represented as a python class, and it
    receives a channel identifier (ID). This test calls the respective object
    for the channel provided by the user (1, 2, 3, or 4) and executs its
    function get_channel_ID(). The returned results should be in agreement with
    the provided channel.

    - Create Image Name: tess of the function that creates the image name based
    on the provided operation mode of the CCD

    - Configura Gain: tess of the function that sets the CCD gain based
    on the provided operation mode of the CCD


Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

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


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        star_magnitude=15,
        sky_magnitude=10,
        ccd_operation_mode=dic,
        channel=1,
        gaussian_std=3,
        star_coordinates=(100, 100),
        bias_level=500,
        sparc4_operation_mode="phot",
        image_dir="a",
    )


# -------------------- Testing the AIS class --------------------------------


def test_star_magnitude_positive_value(ais):
    assert ais.star_magnitude == 15


def test_sky_magnitude_positive_value(ais):
    assert ais.sky_magnitude == 10


def test_channel_value(ais):
    assert ais.channel == 1


def test_gaussian_std_positive_value(ais):
    assert ais.gaussian_std == 3


def test_star_coordinates_value(ais):
    assert ais.star_coordinates == (100, 100)


def test_bias_level(ais):
    assert ais.bias_level == 500


def test_sparc4_operation_mode(ais):
    assert ais.sparc4_operation_mode == "phot"


def test_image_dir(ais):
    assert ais.image_dir == "a\\"


def test_CHC(ais):
    var = 0
    if ais.CHC:
        var = 1
    assert var == 1


def test_SC(ais):
    var = 0
    if ais.SC:
        var = 1
    assert var == 1


def test_TSR(ais):
    var = 0
    if ais.TSR:
        var = 1
    assert var == 1


def test_ASR(ais):
    var = 0
    if ais.ASR:
        var = 1
    assert var == 1


# -----------------Provide a string value to the parameters--------------------


def test_star_magnitude_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator("a", 10, dic)


def test_sky_magnitude_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(20, "a", dic)


def test_gaussian_stddev_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(20, 10, dic, gaussian_std="a")


def test_star_coordinates_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(20, 10, dic, star_coordinates=("a", 10))


def test_bias_level_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, bias_level="a")


def test_sparc4_operation_mode_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, image_dir=1)


def test_image_dir_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, image_dir=1)


# ------------------provide a negative value to the parameters-----------------


def test_star_magnitude_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(-100, 10, dic)


def test_sky_magnitude_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, -10, dic)


def test_gaussian_stddev_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, gaussian_std=-1)


def test_star_coordinates_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, star_coordinates=(-1, -1))


def test_bias_level_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, bias_level=-1)


# ----------------provide zero to the parameters----------------------


def test_star_magnitude_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(0, 10, dic)


def test_sky_magnitude_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 0, dic)


def test_gaussian_stddev_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, gaussian_std=0)


def test_star_coordinates_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic, star_coordinates=(0, 0))


def test_bias_level_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, bias_level=0)


# -------  provide a wrong parameter to the CCD operation mode---------------


def test_ccd_operation_mode_wrong_keyword_em_mode():
    dic = {
        "em_modee": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_em_gain():
    dic = {
        "em_mode": 0,
        "em_gainn": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_preamp():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preampp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_hss():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hsss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_bin():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "bin": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_t_exp():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_expp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_ccd_temp():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_tempp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyword_image_size():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_sizee": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


# ------------ provide a wrong value to the CCD operation mode-------------


def test_ccd_operation_mode_wrong_keyvalue_em_mode():
    dic = {
        "em_mode": 2,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_1():
    dic = {
        "em_mode": 0,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_2():
    dic = {
        "em_mode": 1,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_3():
    dic = {
        "em_mode": 1,
        "em_gain": "a",
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_preamp():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 3,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_hss():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 0,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_bin():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 3,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_1():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 0,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_2():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_ccd_temp_1():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 0,
        "ccd_temp": 50,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_ccd_temp_2():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_1():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": -1,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_2():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 0,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_3():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 1.1,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


# ---------------------------miscelaneous-------------------------------------


def test_gaussian_stddev_float_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(-100, 10, dic, gaussian_std=3.0)


def test_channel_wrong_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(20, 10, dic, channel=5)


# ---------------------------Missing CCD operation mode parameter -------------


def test_ccd_operation_mode_missing_parameter_em_mode():
    dic = {
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_em_gain():
    dic = {
        "em_mode": 0,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_preamp():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_hss():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_t_exp():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_bin():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_ccd_temp():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


def test_ccd_operation_mode_missing_parameter_image_size():
    dic = {
        "em_mode": 0,
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, dic)


# ----------------------- Channels ID --------------------------------


def test_Channel_output_1():
    ais = Artificial_Image_Simulator(100.0, 10.0, dic, channel=1)
    assert ais.get_channel_ID() == "Channel 1"


def test_Channel_output_2():
    ais = Artificial_Image_Simulator(100.0, 10.0, dic, channel=2)
    assert ais.get_channel_ID() == "Channel 2"


def test_Channel_output_3():
    ais = Artificial_Image_Simulator(100.0, 10.0, dic, channel=3)
    assert ais.get_channel_ID() == "Channel 3"


def test_Channel_output_4():
    ais = Artificial_Image_Simulator(100.0, 10.0, dic, channel=4)
    assert ais.get_channel_ID() == "Channel 4"


# ----------------------- Calculate dark current  -------------------------


def test_calculate_dark_current(ais):
    ais._calculate_dark_current()
    assert round(ais.dark_current, 7) == 5.86e-5


# -------------------------Calculate Read Noise -------------------------


def test_calculate_read_noise(ais):
    ais._calculate_read_noise(dic)
    assert ais.read_noise == 6.67


# ----------------------- Apply spectruns ---------------------------


def test_calcualte_star_spectrum(ais):
    ais._calculate_star_spectrum()
    assert ais.star_spectrum == [100]


def test_calculate_sky_spectrum(ais):
    ais._calculate_sky_spectrum()
    assert ais.sky_spectrum == [100]


def test_apply_atmosphere_spectral_response(ais):
    ais.apply_atmosphere_spectral_response()
    assert ais.star_spectrum == [100]


def test_apply_telescope_spectral_response(ais):
    ais.apply_atmosphere_spectral_response()
    assert ais.star_spectrum == [100]


def test_apply_sparc4_spectral_response(ais):
    ais.apply_sparc4_spectral_response()
    assert ais.star_spectrum == [100]


# -----------------------------test _create_image_name------------------------


@pytest.mark.parametrize(
    "em_mode, em_gain, hss, preamp, binn, t_exp, include_star_mag, image_name",
    [
        (0, 1, 0.1, 1, 1, 1, False, "CONV_HSS0.1_PA1_B1_TEXP1_G1"),
        (0, 1, 0.1, 1, 2, 1, False, "CONV_HSS0.1_PA1_B2_TEXP1_G1"),
        (0, 1, 0.1, 2, 1, 1, False, "CONV_HSS0.1_PA2_B1_TEXP1_G1"),
        (0, 1, 0.1, 2, 2, 1, False, "CONV_HSS0.1_PA2_B2_TEXP1_G1"),
        (0, 1, 1, 1, 1, 1, False, "CONV_HSS1_PA1_B1_TEXP1_G1"),
        (0, 1, 1, 1, 2, 1, False, "CONV_HSS1_PA1_B2_TEXP1_G1"),
        (0, 1, 1, 2, 1, 1, False, "CONV_HSS1_PA2_B1_TEXP1_G1"),
        (0, 1, 1, 2, 2, 1, False, "CONV_HSS1_PA2_B2_TEXP1_G1"),
        (1, 2, 1, 1, 1, 1, False, "EM_HSS1_PA1_B1_TEXP1_G2"),
        (1, 2, 1, 1, 2, 1, False, "EM_HSS1_PA1_B2_TEXP1_G2"),
        (1, 2, 1, 2, 1, 1, False, "EM_HSS1_PA2_B1_TEXP1_G2"),
        (1, 2, 1, 2, 2, 1, False, "EM_HSS1_PA2_B2_TEXP1_G2"),
        (1, 2, 10, 1, 1, 1, False, "EM_HSS10_PA1_B1_TEXP1_G2"),
        (1, 2, 10, 1, 2, 1, False, "EM_HSS10_PA1_B2_TEXP1_G2"),
        (1, 2, 10, 2, 1, 1, False, "EM_HSS10_PA2_B1_TEXP1_G2"),
        (1, 2, 10, 2, 2, 1, False, "EM_HSS10_PA2_B2_TEXP1_G2"),
        (1, 2, 20, 1, 1, 1, False, "EM_HSS20_PA1_B1_TEXP1_G2"),
        (1, 2, 20, 1, 2, 1, False, "EM_HSS20_PA1_B2_TEXP1_G2"),
        (1, 2, 20, 2, 1, 1, False, "EM_HSS20_PA2_B1_TEXP1_G2"),
        (1, 2, 20, 2, 2, 1, False, "EM_HSS20_PA2_B2_TEXP1_G2"),
        (1, 2, 30, 1, 1, 1, False, "EM_HSS30_PA1_B1_TEXP1_G2"),
        (1, 2, 30, 1, 2, 1, False, "EM_HSS30_PA1_B2_TEXP1_G2"),
        (1, 2, 30, 2, 1, 1, False, "EM_HSS30_PA2_B1_TEXP1_G2"),
        (1, 2, 30, 2, 2, 1, False, "EM_HSS30_PA2_B2_TEXP1_G2"),
        (1, 2, 30, 2, 2, 2, False, "EM_HSS30_PA2_B2_TEXP2_G2"),
        (1, 2, 30, 2, 2, 1, True, "EM_HSS30_PA2_B2_TEXP1_G2_S15"),
    ],
)
def test_create_image_name(
    ais, em_mode, em_gain, hss, preamp, binn, t_exp, include_star_mag, image_name
):
    dic = {
        "em_mode": em_mode,
        "em_gain": em_gain,
        "preamp": preamp,
        "hss": hss,
        "binn": binn,
        "t_exp": t_exp,
    }
    ais._configure_image_name(ccd_operation_mode=dic, include_star_mag=include_star_mag)
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
    ais._configure_gain(dic)
    assert ais.ccd_gain == ccd_gain


# --------------------------- test create artificial image ----------------


def test_create_artificial_image():
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
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
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
        "image_size": 200,
    }
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
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
        "image_size": 200,
    }
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
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
        "image_size": 200,
    }
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
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
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
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
    ais = Artificial_Image_Simulator(
        100, 10, dic, image_dir=r"C:\Users\denis\Desktop\FITS"
    )
    ais.create_flat_image()
