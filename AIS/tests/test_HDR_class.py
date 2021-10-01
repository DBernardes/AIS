# -*- coding: utf-8 -*-
"""Header class tests.

This script tests the operation of the Header Class. 

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import os

import numpy as np
import pandas as pd
import pytest
from AIS.Header import Header

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


@pytest.fixture
def hdr():
    hdr = Header(dic, 3.0, 9914)
    hdr.create_header()
    return hdr


file = os.path.join("Header", "header.csv")
ss = pd.read_csv(file)
keywords = ss["Keyword"]
comments = ss["Comment"]


# -----------------------teste cria classe -----------------------------


def test_em_mode(hdr):
    assert hdr.em_mode == 0


def test_noise_factor_1(hdr):
    assert hdr.noise_factor == 1


def test_noise_factor_2():
    dic = {
        "em_mode": 1,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 1024,
    }
    hdr = Header(dic, 3, 9914)
    assert hdr.noise_factor == 1.41


def test_em_gain_1(hdr):
    assert hdr.em_gain == 1


def test_em_gain_2(hdr):
    dic = {
        "em_mode": 1,
        "em_gain": 2,
        "preamp": 1,
        "hss": 1,
        "binn": 1,
        "t_exp": 1,
        "ccd_temp": -70,
        "image_size": 1024,
    }
    hdr = Header(dic, 3, 9914)
    assert hdr.em_gain == 2


def test_preamp(hdr):
    assert hdr.preamp == 1


def test_hss(hdr):
    assert hdr.hss == 1


def test_bin(hdr):
    assert hdr.binn == 1


def test_t_exp(hdr):
    assert hdr.t_exp == 1


def test_ccd_gain(hdr):
    assert hdr.ccd_gain == 3.0


def test_serial_number(hdr):
    assert hdr.serial_number == 9914


# ---------------------------- teste read spreadsheet ----------------------


def test_read_spreadsheet(hdr):
    hdr._read_spreadsheet()
    assert hdr.keywords.all() == keywords.all()
    assert hdr.comments.all() == comments.all()


# -----------------------------test _create_image_header---------------------


def test_NAXIS1(hdr):
    assert hdr.hdr["NAXIS1"] == dic["image_size"]


def test_NAXIS2(hdr):
    assert hdr.hdr["NAXIS2"] == dic["image_size"]


def test_HBIN(hdr):
    assert hdr.hdr["HBIN"] == dic["binn"]


def test_VBIN1(hdr):
    assert hdr.hdr["VBIN"] == dic["binn"]


def test_EXPOSURE(hdr):
    assert hdr.hdr["EXPOSURE"] == dic["t_exp"]


def test_READTIME(hdr):
    assert hdr.hdr["READTIME"] == f"{float(dic['hss'])}E-006"


def test_TEMP(hdr):
    assert hdr.hdr["TEMP"] == dic["ccd_temp"]


def test_GAIN(hdr):
    assert hdr.hdr["GAIN"] == 3.0


def test_OUTPTAMP(hdr):
    hdr.create_header()
    assert hdr.hdr["OUTPTAMP"] == "Conventional"


def test_EM_GAIN(hdr):
    assert hdr.hdr["EMGAIN"] == dic["em_gain"]


def test_PREAMP(hdr):
    assert hdr.hdr["PREAMP"] == f"{dic['preamp']}x"


def test_SERNO(hdr):
    hdr.create_header()
    assert hdr.hdr["SERNO"] == 9914
