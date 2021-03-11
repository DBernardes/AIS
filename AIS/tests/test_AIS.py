# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 13:37:22 2021

@author: denis
"""

# -------------------- Testing the AIS class --------------------------------
from AIS import Artificial_Images_Simulator
import pytest


@pytest.fixture
def ais():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    return Artificial_Images_Simulator(100.0, 10.0, 3, dic)


dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1, 'hss': 1, 'bin': 1, 't_exp': 1}


def test_star_flux_positive_value(ais):
    assert ais.star_flux == 100


def test_sky_flux_positive_value(ais):
    assert ais.sky_flux == 10


def test_gaussian_stddev_positive_value(ais):
    assert ais.gaussian_stddev == 3


def test_ccd_temperature():
    ais = Artificial_Images_Simulator(100, 10, 3, dic, ccd_temp=-70)
    assert ais.ccd_temp == -70


def test_bias_level():
    ais = Artificial_Images_Simulator(100, 10, 3, dic, bias_level=500)
    assert ais.bias_level == 500


def test_image_dir():
    ais = Artificial_Images_Simulator(100, 10, 3, dic, image_dir='')
    assert ais.image_dir == ''

# --------------------------------------------------------------------------


def test_star_flux_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator('a', 10, 3, dic)


def test_sky_flux_isnot_a_number():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 'a', 3, dic)


def test_gaussian_stddev_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(-100, 10, 'a', dic)


def test_ccd_temperature_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 3, dic, ccd_temp='a')


def test_bias_level_isnot_number():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 3, dic, bias_level='a')


def test_image_dir_isnot_string():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 3, dic, image_dir=1)
# --------------------------------------------------------------------------


def test_gaussian_stddev_float_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(-100, 10, 3.0, dic)


def test_operation_mode(ais):
    assert ais.ccd_operation_mode == {
        'em_mode': 0, 'em_gain': 2, 'preamp': 1,
        'hss': 1, 'bin': 1, 't_exp': 1}


def test_ccd_temperature_out_of_range():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 3, dic, ccd_temp=-80)


def test_ccd_operation_mode_missing_parameter():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)

# --------------------------------------------------------------------------


def test_star_flux_negative_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(-100, 10, 3, dic)


def test_sky_flux_negative_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, -10, 3, dic)


def test_gaussian_stddev_negative_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, -3, dic)


def test_bias_level_negative_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, -3, dic, bias_level=-1)


# --------------------------------------------------------------------------


def test_star_flux_zero_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(0, 10, 3, dic)


def test_sky_flux_zero_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 0, 3, dic)


def test_gaussian_stddev_zero_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_bias_level_zero_value():
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic, bias_level=0)

# --------------------------------------------------------------------------


def test_ccd_operation_mode_wrong_keyword_1():
    dic = {'em_modee': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyword_2():
    dic = {'em_mode': 0, 'em_gainn': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyword_3():
    dic = {'em_mode': 0, 'em_gain': 2, 'preampp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyword_4():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hsss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyword_5():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'binn': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyword_6():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_expp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)

# --------------------------------------------------------------------------


def test_ccd_operation_mode_wrong_keyvalue_1():
    dic = {'em_mode': 2, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyvalue_2():
    dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyvalue_3():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 3,
           'hss': 1, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyvalue_4():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 0, 'bin': 1, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyvalue_5():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 3, 't_exp': 1}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)


def test_ccd_operation_mode_wrong_keyvalue_6():
    dic = {'em_mode': 0, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 0}
    with pytest.raises(ValueError):
        Artificial_Images_Simulator(100, 10, 0, dic)
