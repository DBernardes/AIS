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
    CCD are provided as a python ccd_operation_modetionary. It is tested if the keyword values
    of the provided ccd_operation_modetionary are correct; it is tested if the values of the
    parameters are in agreement with the manufacture manual.

    - Channel ID: each SPARC4 channel is represented as a python class, and it
    receives a channel identifier (ID). This test calls the respective object
    for the channel provided by the user (1, 2, 3, or 4) and executs its
    function get_channel_ID(). The returned results should be in agreement with
    the provided channel. 


Created on Fri Apr 16 09:10:51 2021

@author: denis
"""


import pytest
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

from .AIS_spectral_response_curves import (
    air_mass,
    ccd_operation_mode,
    l_final,
    l_init,
    l_step,
    magnitude,
    sky_condition,
    sparc4_operation_mode,
    star_temperature,
)

channel = 1
star_coordinates = (100, 100)
image_dir = "a"
bias_level = 500
seeing = 1.5


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=ccd_operation_mode,
        channel=channel,
        star_coordinates=star_coordinates,
        bias_level=bias_level,
        sparc4_operation_mode=sparc4_operation_mode,
        image_dir=image_dir,
        wavelength_interval=(l_init, l_final, l_step),
        star_temperature=star_temperature,
        star_magnitude=magnitude,
        seeing=seeing,
        air_mass=air_mass,
        sky_condition=sky_condition,
    )


# -------------------- Testing the AIS class --------------------------------


def test_channel_value(ais):
    assert ais.channel == channel


def test_star_coordinates_value(ais):
    assert ais.star_coordinates == star_coordinates


def test_bias_level(ais):
    assert ais.bias_level == bias_level


def test_image_dir(ais):
    assert ais.image_dir == image_dir


def test_wavelength_interval(ais):
    assert ais.wavelength_interval == range(400, 1150, 50)


def test_star_temperature(ais):
    assert ais.star_temperature == star_temperature


def test_star_magnitude(ais):
    assert ais.star_magnitude == magnitude


def test_seeing(ais):
    assert ais.seeing == seeing


def test_air_mass(ais):
    assert ais.air_mass == air_mass


def test_sky_condition(ais):
    assert ais.sky_condition == sky_condition


def test_CHC(ais):
    var = 0
    if ais.chc:
        var = 1
    assert var == 1


def test_SC(ais):
    var = 0
    if ais.sc:
        var = 1
    assert var == 1


def test_TSR(ais):
    var = 0
    if ais.tsr:
        var = 1
    assert var == 1


def test_ASR(ais):
    var = 0
    if ais.asr:
        var = 1
    assert var == 1


# -----------------Provide a string value to the parameters--------------------


def test_star_coordinates_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_coordinates=("a", 10))


def test_bias_level_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, bias_level="a")


def test_image_dir_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, image_dir=1)


def test_wavelength_interval_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, wavelength_interval="a")


def test_star_temperature_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_temperature="a")


def test_star_magnitude_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_magnitude="a")


def test_seeing_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, seeing="a")


def test_air_mass_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, air_mass="a")


# ------------------provide a negative value to the parameters-----------------


def test_star_coordinates_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_coordinates=(-1, -1))


def test_bias_level_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, bias_level=-1)


def test_wavelength_interval_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, wavelength_interval=(-20, 20, 20)
        )


def test_star_magnitude_negative_value():
    ais = Artificial_Image_Simulator(ccd_operation_mode, star_magnitude=-1)
    assert ais.star_magnitude == -1


def test_seeing_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, seeing=-1.5)


def test_air_mass_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, air_mass=-air_mass)


# ----------------provide zero to the parameters----------------------


def test_star_coordinates_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_coordinates=(0, 0))


def test_bias_level_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, bias_level=0)


def test_wavelength_interval_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, wavelength_interval=(0, 20, 20))


def test_star_temperature_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_temperature=0)


def test_star_magnitude_zero_value():
    ais = Artificial_Image_Simulator(ccd_operation_mode, star_magnitude=0)
    assert ais.star_magnitude == 0


def test_seeing_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, seeing=0)


def test_air_mass_zero_value():
    ais = Artificial_Image_Simulator(ccd_operation_mode, air_mass=0)
    assert ais.air_mass == 0


def test_sky_condition_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, sky_condition=0)


# -------  provide a wrong parameter to the CCD operation mode---------------


def test_ccd_operation_mode_wrong_keyword_em_mode():
    ccd_operation_mode = {
        "em_modee": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_em_gain():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gainn": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_preamp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preampp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_hss():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hsss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_bin():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "bin": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_t_exp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_expp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_ccd_temp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_tempp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyword_image_size():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_sizee": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


# ------------ provide a wrong value to the CCD operation mode-------------


def test_ccd_operation_mode_wrong_keyvalue_em_mode():
    ccd_operation_mode = {
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
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_1():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_2():
    ccd_operation_mode = {
        "em_mode": "EM",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_3():
    ccd_operation_mode = {
        "em_mode": "EM",
        "em_gain": "a",
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_preamp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 3,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_hss():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 0,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_bin():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 3,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_1():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 0,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_2():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_ccd_temp_1():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 0,
        "ccd_temp": 50,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_ccd_temp_2():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_1():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": -1,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_2():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 0,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_wrong_keyvalue_t_image_size_3():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": "a",
        "ccd_temp": -70,
        "image_size": 1.1,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


# -----------------------Test SPARC4 operation mode ---------------------------------------


def test_dic_without_acquisition_mode():
    sparc4_operation_mode = {}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            sparc4_operation_mode=sparc4_operation_mode,
            ccd_operation_mode=ccd_operation_mode,
        )


def test_wrong_acquisition_mode():
    dic = {"acquisition_mode": "a"}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            sparc4_operation_mode=dic, ccd_operation_mode=ccd_operation_mode
        )


def test_acquisition_photometric(ais):
    assert ais.sparc4_operation_mode["acquisition_mode"] == "photometric"


def test_acquisition_polarimetric(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["acquisition_mode"] == "polarimetric"


def test_acquisition_photometric_with_unnecessary_keywords(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "photometric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
        )


def test_acquisition_polarimetric_and_wrong_keyword(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
        "blah": 1,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
        )


def test_acquisition_polarimetric_missing_keyword(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
        )


def test_acquisition_polarimetric_calibration_wheel_polarizer(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "polarizer",
        "retarder": "half",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["calibration_wheel"] == "polarizer"


def test_acquisition_polarimetric_calibration_wheel_depolarizer(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "depolarizer",
        "retarder": "half",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["calibration_wheel"] == "depolarizer"


def test_acquisition_polarimetric_calibration_wheel_empty(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "half",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["calibration_wheel"] == "empty"


def test_acquisition_polarimetric_calibration_wheel_wrong_value(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "blah",
        "retarder": "half",
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
        )


def test_acquisition_polarimetric_retarder_qaurter(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "polarizer",
        "retarder": "quarter",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["retarder"] == "quarter"


def test_acquisition_polarimetric_retarder_half(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "depolarizer",
        "retarder": "half",
    }
    ais = Artificial_Image_Simulator(
        ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
    )
    assert ais.sparc4_operation_mode["retarder"] == "half"


def test_acquisition_polarimetric_retarder_wrong_value(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "polarizer",
        "retarder": "blah",
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(
            ccd_operation_mode, sparc4_operation_mode=sparc4_operation_mode
        )


# ---------------------------miscelaneous--------------------------------------------------


def test_channel_wrong_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, channel=5)


def test_star_coordinates_greater_than_image_size():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode, star_coordinates=(2000, 2000))


# ---------------------------Missing CCD operation mode parameter -------------


def test_ccd_operation_mode_missing_parameter_em_mode():
    ccd_operation_mode = {
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_em_gain():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_preamp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_hss():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_t_exp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_bin():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_ccd_temp():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "image_size": 200,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


def test_ccd_operation_mode_missing_parameter_image_size():
    ccd_operation_mode = {
        "em_mode": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
    }
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(ccd_operation_mode)


# ----------------------- Channels ID --------------------------------


def test_Channel_output_1():
    ais = Artificial_Image_Simulator(ccd_operation_mode, channel=1)
    assert ais.get_channel_id() == "Channel 1"


def test_Channel_output_2():
    ais = Artificial_Image_Simulator(ccd_operation_mode, channel=2)
    assert ais.get_channel_id() == "Channel 2"


def test_Channel_output_3():
    ais = Artificial_Image_Simulator(ccd_operation_mode, channel=3)
    assert ais.get_channel_id() == "Channel 3"


def test_Channel_output_4():
    ais = Artificial_Image_Simulator(ccd_operation_mode, channel=4)
    assert ais.get_channel_id() == "Channel 4"
