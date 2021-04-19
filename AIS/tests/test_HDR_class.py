# -*- coding: utf-8 -*-
"""HDR class tests.

This script tests the operation of the Header Class. 

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


from HDR import Header
import pytest

dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
       'hss': 1, 'bin': 1, 't_exp': 1, 'ccd_temp': -70}


@pytest.fixture
def hdr():
    return Header(dic, 3.0, 9914)


# -----------------------teste cria classe -----------------------------


def test_em_mode(hdr):
    assert hdr.em_mode == 0


def test_noise_factor_1(hdr):
    assert hdr.noise_factor == 1


def test_noise_factor_2():
    dic = {'em_mode': 1, 'em_gain': 2, 'preamp': 1,
           'hss': 1, 'bin': 1, 't_exp': 1, 'ccd_temp': -70}
    hdr = Header(dic, 3, 9914)
    assert hdr.noise_factor == 1.41


def test_em_gain_1(hdr):
    assert hdr.em_gain == 1


def test_em_gain_2(hdr):
    hdr.em_mode = 1
    hdr.em_gain = 2
    assert hdr.em_gain == 2


def test_preamp(hdr):
    assert hdr.preamp == 1


def test_hss(hdr):
    assert hdr.hss == 1


def test_bin(hdr):
    assert hdr.bin == 1


def test_t_exp(hdr):
    assert hdr.t_exp == 1


def test_ccd_gain(hdr):
    assert hdr.ccd_gain == 3.0


def test_serial_number(hdr):
    assert hdr.serial_number == 9914

# -----------------------------test _create_image_header---------------------


def test_NAXIS1(hdr):
    hdr._create_header()
    assert hdr.hdr['NAXIS1'] == 200


def test_NAXIS2(hdr):
    hdr._create_header()
    assert hdr.hdr['NAXIS2'] == 200


def test_HBIN(hdr):
    hdr._create_header()
    assert hdr.hdr['HBIN'] == 1


def test_VBIN1(hdr):
    hdr._create_header()
    assert hdr.hdr['VBIN'] == 1


def test_EXPOSURE(hdr):
    hdr._create_header()
    assert hdr.hdr['EXPOSURE'] == 1


def test_TEMP(hdr):
    hdr._create_header()
    assert hdr.hdr['READTIME'] == '1.0E-006'


def test_GAIN(hdr):
    hdr._create_header()
    assert hdr.hdr['GAIN'] == 3.0


def test_OUTPTAMP(hdr):
    hdr._create_header()
    assert hdr.hdr['OUTPTAMP'] == 'Conventional'


def test_SERNO(hdr):
    hdr._create_header()
    assert hdr.hdr['SERNO'] == 9914
