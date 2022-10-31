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
    moon_condition,
    sky_condition,
    sparc4_operation_mode,
)


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(
        ccd_operation_mode=ccd_operation_mode)


# # -------------------- Testing the AIS class --------------------------------

def test_ccd_operation_mode(ais):
    assert ais.ccd_operation_mode == ccd_operation_mode

# -------  provide a wrong parameter to the CCD operation mode---------------


def test_ccd_operation_mode_wrong_keyword_em_mode():
    ccd_operation_mode = {
        "em_modee": "Conv",
        "em_gain": 1,
        "preamp": 1,
        "readout": 1,
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


def test_dic_without_acquisition_mode(ais):
    sparc4_operation_mode = {}
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode,
        )


def test_wrong_acquisition_mode(ais):
    dic = {"acquisition_mode": "a"}
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(sparc4_operation_mode=dic)


def test_acquisition_photometric_with_unnecessary_keywords(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "photometric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
    }
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode)


def test_acquisition_polarimetric_and_wrong_keyword(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
        "retarder": "quarter",
        "blah": 1,
    }
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode)


def test_acquisition_polarimetric_missing_keyword(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "empty",
    }
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode)


def test_acquisition_polarimetric_calibration_wheel_wrong_value(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "blah",
        "retarder": "half",
    }
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode)


def test_acquisition_polarimetric_retarder_wrong_value(ais):
    sparc4_operation_mode = {
        "acquisition_mode": "polarimetric",
        "calibration_wheel": "polarizer",
        "retarder": "blah",
    }
    with pytest.raises(ValueError):
        ais.apply_sparc4_spectral_response(
            sparc4_operation_mode=sparc4_operation_mode)


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
