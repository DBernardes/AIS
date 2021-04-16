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

# -----------------provide a string value to the parameters--------------------


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

# ---------------------------miscelaneous-------------------------------------


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

# ------------------provide a negative value to the parameters-----------------


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


# ----------------provide zero to the parameters----------------------


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

# -------  provide a wrong parameter to the CCD operation mode---------------


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

# ------------ provide a wrong value to the CCD operation mode-------------


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

# -----------------------função _write_image_mode-----------------------------


def test_write_image_mode_em_mode(ais):
    ais._write_image_mode()
    assert ais.em_mode == 0


def test_write_image_mode_noise_factor_1(ais):
    ais._write_image_mode()
    assert ais.noise_factor == 1


def test_write_image_mode_noise_factor_2(ais):
    dic = {'em_mode': 1, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    ais.ccd_operation_mode = dic
    ais._write_image_mode()
    assert ais.noise_factor == 1.41


def test_write_image_mode_em_gain_1(ais):
    ais._write_image_mode()
    assert ais.em_gain == 1


def test_write_image_mode_em_gain_2(ais):
    dic = {'em_mode': 1, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1}
    ais.ccd_operation_mode = dic
    ais._write_image_mode()
    assert ais.em_gain == 2


def test_write_image_mode_preamp(ais):
    ais._write_image_mode()
    assert ais.noise_factor == 1


def test_write_image_mode_hss(ais):
    ais._write_image_mode()
    assert ais.noise_factor == 1


def test_write_image_mode_bin(ais):
    ais._write_image_mode()
    assert ais.noise_factor == 1


def test_write_image_mode_t_exp(ais):
    ais._write_image_mode()
    assert ais.noise_factor == 1

# ---------------------_calc_dark_current--------------------------------


def test_dark_current_9914(ais):
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(5.85975599e-5, 7)


def test_dark_current_9915(ais):
    ais.serial_number = 9915
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(1.466809375e-4, 7)


def test_dark_current_9916(ais):
    ais.serial_number = 9916
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(8.688095666e-5, 7)


def test_dark_current_9917(ais):
    ais.serial_number = 9917
    ais._calc_dark_current()
    assert round(ais.dark_current, 7) == round(2.313304035e-4, 7)

# -----------------------------_calc_read_noise--------------------------


@pytest.mark.parametrize(
    'em_mode, em_gain, hss, preamp, binn, read_noise',
    [(0, 2, 0.1, 1, 1, 8.78),
     (0, 2, 0.1, 1, 2, 8.84),
     (0, 2, 0.1, 2, 1, 3.46),
     (0, 2, 0.1, 2, 2, 3.27),
     (0, 2, 1, 1, 1, 6.67),
     (0, 2, 1, 1, 2, 6.94),
     (0, 2, 1, 2, 1, 4.76),
     (0, 2, 1, 2, 2, 4.79),
     (1, 2, 1, 1, 1, 24.64),
     (1, 2, 1, 1, 2, 33.76),
     (1, 2, 1, 2, 1, 12.05),
     (1, 2, 1, 2, 2, 14.55),
     (1, 2, 10, 1, 1, 83.68),
     (1, 2, 10, 1, 2, 82.93),
     (1, 2, 10, 2, 1, 41.71),
     (1, 2, 10, 2, 2, 41.82),
     (1, 2, 20, 1, 1, 160.06),
     (1, 2, 20, 1, 2, 161.98),
     (1, 2, 20, 2, 1, 66.01),
     (1, 2, 20, 2, 2, 72.71),
     (1, 2, 30, 1, 1, 262.01),
     (1, 2, 30, 1, 2, 273.19),
     (1, 2, 30, 2, 1, 169.25),
     (1, 2, 30, 2, 2, 143.59),
     ]
)
def test_calc_read_noise(ais, em_mode, em_gain, hss, preamp, binn, read_noise):
    ais.em_mode = em_mode
    ais.em_gain = em_gain
    ais.hss = hss
    ais.preamp = preamp
    ais.preamp = preamp
    ais.bin = binn
    ais._calc_read_noise()
    assert round(ais.read_noise, 2) == read_noise

# -----------------------------test _create_image_name------------------------


@pytest.mark.parametrize(
    'em_mode, em_gain, hss, preamp, binn, t_exp, include_star_flux, image_name',
    [(0, 2, 0.1, 1, 1, 1, False, 'CONV_HSS0.1_PA1_B1_TEXP1_G1'),
     (0, 2, 0.1, 1, 2, 1, False, 'CONV_HSS0.1_PA1_B2_TEXP1_G1'),
     (0, 2, 0.1, 2, 1, 1, False, 'CONV_HSS0.1_PA2_B1_TEXP1_G1'),
     (0, 2, 0.1, 2, 2, 1, False, 'CONV_HSS0.1_PA2_B2_TEXP1_G1'),
     (0, 2, 1, 1, 1, 1, False, 'CONV_HSS1_PA1_B1_TEXP1_G1'),
     (0, 2, 1, 1, 2, 1, False, 'CONV_HSS1_PA1_B2_TEXP1_G1'),
     (0, 2, 1, 2, 1, 1, False, 'CONV_HSS1_PA2_B1_TEXP1_G1'),
     (0, 2, 1, 2, 2, 1, False, 'CONV_HSS1_PA2_B2_TEXP1_G1'),
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
        ais, em_mode, em_gain, hss, preamp, binn, t_exp, include_star_flux, image_name):
    dic = {'em_mode': em_mode, 'em_gain': em_gain, 'preamp': preamp,
           'hss': hss, 'bin': binn, 't_exp': t_exp}
    ais.ccd_operation_mode = dic
    ais._write_image_mode()
    ais._configure_image_name(include_star_flux=include_star_flux)
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
    ais.em_mode = em_mode
    ais.em_gain = em_gain
    ais.hss = hss
    ais.preamp = preamp
    ais.preamp = preamp
    ais.bin = binn
    ais._configure_gain()
    assert ais.gain == ccd_gain


# -----------------------------test _create_image_header---------------------

def test_NAXIS1(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['NAXIS1'] == 200


def test_NAXIS2(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['NAXIS2'] == 200


def test_HBIN(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['HBIN'] == 1


def test_VBIN1(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['VBIN'] == 1


def test_EXPOSURE(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['EXPOSURE'] == 1


def test_TEMP(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['READTIME'] == '1.0E-006'


def test_GAIN(ais):
    ais._write_image_mode()
    ais._configure_gain()
    ais._create_image_header()
    assert ais.hdr['GAIN'] == 3.37


def test_OUTPTAMP(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['OUTPTAMP'] == 'Conventional'


def test_SERNO(ais):
    ais._write_image_mode()
    ais._create_image_header()
    assert ais.hdr['SERNO'] == 9914
