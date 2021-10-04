# -*- coding: utf-8 -*-
"""AIS operation tests.

This script presents the tests of the AIS initialization. These tests are:

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


Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


from random import gauss

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

l_init, l_final, l_step = 400, 1100, 50
channel = 1
gaussian_std = 3
star_coordinates = (100, 100)
sparc4_operation_mode = "phot"
image_dir = "a"
star_temperature = 5700
star_magnitude = 22
bias_level = 500


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=dic,
        channel=channel,
        gaussian_std=gaussian_std,
        star_coordinates=star_coordinates,
        bias_level=bias_level,
        sparc4_operation_mode=sparc4_operation_mode,
        image_dir=image_dir,
        wavelength_interval=(l_init, l_final, l_step),
        star_temperature=star_temperature,
        star_magnitude=star_magnitude,
    )


# -------------------- Testing the AIS class --------------------------------


def test_channel_value(ais):
    assert ais.channel == channel


def test_gaussian_std_positive_value(ais):
    assert ais.gaussian_std == gaussian_std


def test_star_coordinates_value(ais):
    assert ais.star_coordinates == star_coordinates


def test_bias_level(ais):
    assert ais.bias_level == bias_level


def test_sparc4_operation_mode(ais):
    assert ais.sparc4_operation_mode == sparc4_operation_mode


def test_image_dir(ais):
    assert ais.image_dir == image_dir


def test_wavelength_interval(ais):
    assert ais.wavelength_interval == range(400, 1150, 50)


def test_star_temperature(ais):
    assert ais.star_temperature == star_temperature


def test_star_magnitude(ais):
    assert ais.star_magnitude == star_magnitude


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


def test_gaussian_stddev_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, gaussian_std="a")


def test_star_coordinates_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_coordinates=("a", 10))


def test_bias_level_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, bias_level="a")


def test_sparc4_operation_mode_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, sparc4_operation_mode=1)


def test_image_dir_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, image_dir=1)


def test_wavelength_interval_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, wavelength_interval="a")


def test_star_temperature_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_temperature="a")


def test_star_magnitude_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_magnitude="a")


# ------------------provide a negative value to the parameters-----------------


def test_gaussian_stddev_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, gaussian_std=-1)


def test_star_coordinates_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_coordinates=(-1, -1))


def test_bias_level_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, bias_level=-1)


def test_wavelength_interval_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, wavelength_interval=(-20, 20, 20))


def test_star_temperature_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_temperature=-1)


def test_star_magnitude_negative_value():
    ais = Artificial_Image_Simulator(dic, star_magnitude=-1)
    assert ais.star_magnitude == -1


# ----------------provide zero to the parameters----------------------


def test_gaussian_stddev_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, gaussian_std=0)


def test_star_coordinates_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_coordinates=(0, 0))


def test_bias_level_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, bias_level=0)


def test_wavelength_interval_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, wavelength_interval=(0, 20, 20))


def test_star_temperature_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, star_temperature=0)


def test_star_magnitude_zero_value():
    ais = Artificial_Image_Simulator(dic, star_magnitude=0)
    assert ais.star_magnitude == 0


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


# ---------------------------miscelaneous-------------------------------------


def test_gaussian_stddev_float_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, gaussian_std=3.0)


def test_channel_wrong_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(dic, channel=5)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


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
        Artificial_Image_Simulator(dic)


# ----------------------- Channels ID --------------------------------


def test_Channel_output_1():
    ais = Artificial_Image_Simulator(dic, channel=1)
    assert ais.get_channel_ID() == "Channel 1"


def test_Channel_output_2():
    ais = Artificial_Image_Simulator(dic, channel=2)
    assert ais.get_channel_ID() == "Channel 2"


def test_Channel_output_3():
    ais = Artificial_Image_Simulator(dic, channel=3)
    assert ais.get_channel_ID() == "Channel 3"


def test_Channel_output_4():
    ais = Artificial_Image_Simulator(dic, channel=4)
    assert ais.get_channel_ID() == "Channel 4"
