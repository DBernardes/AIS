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


from AIS import Artificial_Image_Simulator
import pytest

dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
       'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}


@pytest.fixture
def ais():
    return Artificial_Image_Simulator(star_flux=100.0,
                                      sky_flux=10.0,
                                      gaussian_std=3,
                                      ccd_operation_mode=dic,
                                      channel=1)

# -------------------- Testing the AIS class --------------------------------


def test_star_flux_positive_value(ais):
    assert ais.star_flux == 100


def test_sky_flux_positive_value(ais):
    assert ais.sky_flux == 10


def test_gaussian_std_positive_value(ais):
    assert ais.gaussian_std == 3


def test_bias_level():
    ais = Artificial_Image_Simulator(100, 10, 3, dic, 1, bias_level=300)
    assert ais.bias_level == 300


def test_image_dir():
    ais = Artificial_Image_Simulator(100, 10, 3, dic, 1, image_dir='a')
    assert ais.image_dir == 'a'

# -----------------Provide a string value to the parameters--------------------


def test_star_flux_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator('a', 10, 3, dic, 1)


def test_sky_flux_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 'a', 3, dic, 1)


def test_gaussian_stddev_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(-100, 10, 'a', dic, 1)


def test_bias_level_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 3, dic, 1, bias_level='a')


def test_image_dir_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 3, dic, 1, image_dir=1)


# ------------------provide a negative value to the parameters-----------------


def test_star_flux_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(-100, 10, 3, dic, 1)


def test_sky_flux_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, -10, 3, dic, 1)


def test_gaussian_stddev_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, -3, dic, 1)


def test_bias_level_negative_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, -3, dic, 1, bias_level=-1)

# ----------------provide zero to the parameters----------------------


def test_star_flux_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(0, 10, 3, dic, 1)


def test_sky_flux_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 0, 3, dic, 1)


def test_gaussian_stddev_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_bias_level_zero_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1, bias_level=0)

# -------  provide a wrong parameter to the CCD operation mode---------------


def test_ccd_operation_mode_wrong_keyword_em_mode():
    dic = {'em_modee': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_em_gain():
    dic = {'em_mode': 0, 'em_gainn': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_preamp():
    dic = {'em_mode': 0, 'em_gain': 2, 'preampp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_hss():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hsss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_bin():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_t_exp():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_expp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyword_ccd_temp():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_tempp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)

# ------------ provide a wrong value to the CCD operation mode-------------


def test_ccd_operation_mode_wrong_keyvalue_em_mode():
    dic = {'em_mode': 2, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_1():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_2():
    dic = {'em_mode': 1, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_em_gain_3():
    dic = {'em_mode': 1, 'em_gain': 'a', 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_preamp():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 3,
           'hss': 1, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_hss():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 0, 'binn': 1, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_bin():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 3, 't_exp': 1, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_1():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 0, 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_t_exp_2():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 'a', 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_ccd_temp_1():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 0, 'ccd_temp': 50}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


def test_ccd_operation_mode_wrong_keyvalue_t_ccd_temp_2():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 'a', 'ccd_temp': -70}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


# ---------------------------miscelaneous-------------------------------------


def test_gaussian_stddev_float_value():
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(-100, 10, 3.0, dic, 1)


def test_ccd_operation_mode_missing_parameter():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1}
    with pytest.raises(ValueError):
        Artificial_Image_Simulator(100, 10, 0, dic, 1)


# ----------------------- Channels ID --------------------------------


def test_Channel_output_1():
    ais = Artificial_Image_Simulator(100.0, 10.0, 3, dic, channel=1)
    assert ais.get_channel_ID() == 'Channel 1'


def test_Channel_output_2():
    ais = Artificial_Image_Simulator(100.0, 10.0, 3, dic, channel=2)
    assert ais.get_channel_ID() == 'Channel 2'


def test_Channel_output_3():
    ais = Artificial_Image_Simulator(100.0, 10.0, 3, dic, channel=3)
    assert ais.get_channel_ID() == 'Channel 3'


def test_Channel_output_4():
    ais = Artificial_Image_Simulator(100.0, 10.0, 3, dic, channel=4)
    assert ais.get_channel_ID() == 'Channel 4'

# -----------------------------test _create_image_name------------------------


@pytest.mark.parametrize(
    'em_mode, em_gain, hss, preamp, binn, t_exp, include_star_flux, image_name',
    [(0, 1, 0.1, 1, 1, 1, False, 'CONV_HSS0.1_PA1_B1_TEXP1_G1'),
     (0, 1, 0.1, 1, 2, 1, False, 'CONV_HSS0.1_PA1_B2_TEXP1_G1'),
     (0, 1, 0.1, 2, 1, 1, False, 'CONV_HSS0.1_PA2_B1_TEXP1_G1'),
     (0, 1, 0.1, 2, 2, 1, False, 'CONV_HSS0.1_PA2_B2_TEXP1_G1'),
     (0, 1, 1, 1, 1, 1, False, 'CONV_HSS1_PA1_B1_TEXP1_G1'),
     (0, 1, 1, 1, 2, 1, False, 'CONV_HSS1_PA1_B2_TEXP1_G1'),
     (0, 1, 1, 2, 1, 1, False, 'CONV_HSS1_PA2_B1_TEXP1_G1'),
     (0, 1, 1, 2, 2, 1, False, 'CONV_HSS1_PA2_B2_TEXP1_G1'),
     (1, 2, 1, 1, 1, 1, False, 'EM_HSS1_PA1_B1_TEXP1_G2'),
     (1, 2, 1, 1, 2, 1, False, 'EM_HSS1_PA1_B2_TEXP1_G2'),
     (1, 2, 1, 2, 1, 1, False, 'EM_HSS1_PA2_B1_TEXP1_G2'),
     (1, 2, 1, 2, 2, 1, False, 'EM_HSS1_PA2_B2_TEXP1_G2'),
     (1, 2, 10, 1, 1, 1, False, 'EM_HSS10_PA1_B1_TEXP1_G2'),
     (1, 2, 10, 1, 2, 1, False, 'EM_HSS10_PA1_B2_TEXP1_G2'),
     (1, 2, 10, 2, 1, 1, False, 'EM_HSS10_PA2_B1_TEXP1_G2'),
     (1, 2, 10, 2, 2, 1, False, 'EM_HSS10_PA2_B2_TEXP1_G2'),
     (1, 2, 20, 1, 1, 1, False, 'EM_HSS20_PA1_B1_TEXP1_G2'),
     (1, 2, 20, 1, 2, 1, False, 'EM_HSS20_PA1_B2_TEXP1_G2'),
     (1, 2, 20, 2, 1, 1, False, 'EM_HSS20_PA2_B1_TEXP1_G2'),
     (1, 2, 20, 2, 2, 1, False, 'EM_HSS20_PA2_B2_TEXP1_G2'),
     (1, 2, 30, 1, 1, 1, False, 'EM_HSS30_PA1_B1_TEXP1_G2'),
     (1, 2, 30, 1, 2, 1, False, 'EM_HSS30_PA1_B2_TEXP1_G2'),
     (1, 2, 30, 2, 1, 1, False, 'EM_HSS30_PA2_B1_TEXP1_G2'),
     (1, 2, 30, 2, 2, 1, False, 'EM_HSS30_PA2_B2_TEXP1_G2'),
     (1, 2, 30, 2, 2, 2, False, 'EM_HSS30_PA2_B2_TEXP2_G2'),
     (1, 2, 30, 2, 2, 1, True, 'EM_HSS30_PA2_B2_TEXP1_G2_S100.0'),
     ]
)
def test_create_image_name(
        ais, em_mode, em_gain, hss, preamp, binn,
        t_exp, include_star_flux, image_name):
    dic = {'em_mode': em_mode, 'em_gain': em_gain, 'preamp': preamp,
           'hss': hss, 'binn': binn, 't_exp': t_exp}
    ais._configure_image_name(ccd_operation_mode=dic,
                              include_star_flux=include_star_flux)
    assert ais.image_name == image_name

# -----------------------------test _configure_gain------------------------


@pytest.mark.parametrize(
    'em_mode, em_gain, hss, preamp, binn, ccd_gain',
    [(0, 2, 0.1, 1, 1, 3.35),
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
     ]
)
def test_configure_gain(ais, em_mode, em_gain, hss, preamp, binn, ccd_gain):
    dic = {'em_mode': em_mode, 'em_gain': em_gain, 'preamp': preamp,
           'hss': hss, 'binn': binn, 't_exp': 1, 'ccd_temp': -70}
    ais._configure_gain(dic)
    assert ais.ccd_gain == ccd_gain

# --------------------------- test create artificial image ----------------


def test_create_artificial_image(ais):
    ais.create_artificial_image()
