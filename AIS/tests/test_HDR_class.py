# -*- coding: utf-8 -*-
"""HDR class tests.

This script tests the operation of the Header Class. 

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


from HDR import Header
import pytest

dic = {'em_mode': 0, 'em_gain': 1, 'preamp': 1,
       'hss': 1, 'bin': 1, 't_exp': 1}


@pytest.fixture
def hdr():
    return Header(dic)


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
    hdr._configure_gain()
    hdr._create_header()
    assert hdr.hdr['GAIN'] == 3.37


def test_OUTPTAMP(hdr):
    hdr._create_header()
    assert hdr.hdr['OUTPTAMP'] == 'Conventional'


def test_SERNO(hdr):
    hdr._create_header()
    assert hdr.hdr['SERNO'] == 9914
