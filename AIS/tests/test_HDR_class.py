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

ccd_gain = 3
serial_number = 9914


@pytest.fixture
def hdr():
    hdr = Header(dic, ccd_gain, serial_number)
    hdr.create_header()
    return hdr


file = os.path.join("Header", "header.csv")
ss = pd.read_csv(file)
keywords = ss["Keyword"]
comments = ss["Comment"]


# -----------------------teste cria classe -----------------------------


def test_em_mode(hdr):
    assert hdr.em_mode == dic["em_mode"]


def test_noise_factor_1(hdr):
    assert hdr.NOISE_FACTOR == 1


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
    assert hdr.NOISE_FACTOR == 1.41


def test_em_gain_1(hdr):
    assert hdr.em_gain == dic["em_gain"]


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
    assert hdr.em_gain == dic["em_gain"]


def test_preamp(hdr):
    assert hdr.preamp == dic["preamp"]


def test_hss(hdr):
    assert hdr.hss == dic["hss"]


def test_bin(hdr):
    assert hdr.binn == dic["binn"]


def test_t_exp(hdr):
    assert hdr.t_exp == dic["t_exp"]


def test_ccd_gain(hdr):
    assert hdr.ccd_gain == ccd_gain


def test_serial_number(hdr):
    assert hdr.serial_number == serial_number


# ---------------------------- teste read spreadsheet ----------------------


def test_read_spreadsheet(hdr):
    hdr._read_spreadsheet()
    assert hdr.keywords.all() == keywords.all()
    assert hdr.comments.all() == comments.all()


# -----------------------------test _create_image_header---------------------


def test_BITPIX(hdr):
    assert hdr.hdr["BITPIX"] == "16"


def test_EXTEND(hdr):
    assert hdr.hdr["EXTEND"] == "T"


def test_COMMENT(hdr):
    hdr.hdr["COMMENT"]


def test_FILENAME(hdr):
    assert hdr.hdr["FILENAME"] == "hats-24_I_transito_001"


def test_OBJECT(hdr):
    assert hdr.hdr["OBJECT"] == ""


def test_CHANNEL(hdr):
    assert hdr.hdr["CHANNEL"] == ""


def test_NEXPSEQ(hdr):
    assert hdr.hdr["NEXPSEQ"] == ""


def test_READMODE(hdr):
    assert hdr.hdr["READMODE"] == "Image"


def test_PREAMP(hdr):
    assert hdr.hdr["PREAMP"] == str(dic["preamp"]) + "x"


def test_OBSTYPE(hdr):
    assert hdr.hdr["OBSTYPE"] == ""


def test_VSHIFT(hdr):
    assert hdr.hdr["VSHIFT"] == "4.33E-06"


def test_EMMODE(hdr):
    assert hdr.hdr["EMMODE"] == "Conventional"


def test_VBIN(hdr):
    assert hdr.hdr["VBIN"] == dic["binn"]


def test_INITLIN(hdr):
    assert hdr.hdr["INITLIN"] == ""


def test_KSERLEN(hdr):
    assert hdr.hdr["KSERLEN"] == ""


def test_OPENSHT(hdr):
    assert hdr.hdr["OPENSHT"] == ""


def test_FINALLIN(hdr):
    assert hdr.hdr["FINALLIN"] == ""


def test_DATE(hdr):
    assert hdr.hdr["DATE"] == "2017-07-14T00:00:58"


def test_GAINERR(hdr):
    assert hdr.hdr["GAINERR"] == ""


def test_OBSERVER(hdr):
    assert hdr.hdr["OBSERVER"] == ""


def test_NIGHTDIR(hdr):
    assert hdr.hdr["NIGHTDIR"] == ""


def test_INSTRUME(hdr):
    assert hdr.hdr["INSTRUME"] == ""


def test_SEQINDEX(hdr):
    assert hdr.hdr["SEQINDEX"] == ""


def test_ACQMODE(hdr):
    assert hdr.hdr["ACQMODE"] == "Single"


def test_READOUT(hdr):
    assert hdr.hdr["READOUT"] == str(1 / dic["hss"]) + "E-006"


def test_READMOD(hdr):
    assert hdr.hdr["READMOD"] == ""


def test_TRIGGER(hdr):
    assert hdr.hdr["TRIGGER"] == "External"


def test_EMGAIN(hdr):
    assert hdr.hdr["EMGAIN"] == dic["em_gain"]


def test_HBIN(hdr):
    assert hdr.hdr["HBIN"] == dic["binn"]


def test_INITCOL(hdr):
    assert hdr.hdr["INITCOL"] == ""


def test_SHUTTER(hdr):
    assert hdr.hdr["SHUTTER"] == ""


def test_CLOSESHT(hdr):
    assert hdr.hdr["CLOSESHT"] == ""


def test_FINALCOL(hdr):
    assert hdr.hdr["FINALCOL"] == ""


def test_GAIN(hdr):
    assert hdr.hdr["GAIN"] == ccd_gain


def test_RDNOISE(hdr):
    assert hdr.hdr["RDNOISE"] == ""


def test_PROGRID(hdr):
    assert hdr.hdr["PROGRID"] == ""


def test_RNOISERR(hdr):
    assert hdr.hdr["RNOISERR"] == ""


def test_SKYFIBER(hdr):
    assert hdr.hdr["SKYFIBER"] == ""


def test_CALFIBER(hdr):
    assert hdr.hdr["CALFIBER"] == ""


def test_UTDATE(hdr):
    assert hdr.hdr["UTDATE"] == ""


def test_LOCTIME(hdr):
    assert hdr.hdr["LOCTIME"] == ""


def test_MJD(hdr):
    assert hdr.hdr["MJD"] == ""


def test_OBSLON(hdr):
    assert hdr.hdr["OBSLON"] == ""


def test_OBSLAT(hdr):
    assert hdr.hdr["OBSLAT"] == ""


def test_OBSALT(hdr):
    assert hdr.hdr["OBSALT"] == ""


def test_EXPTIME(hdr):
    assert hdr.hdr["EXPTIME"] == dic["t_exp"]


def test_TELRA(hdr):
    assert hdr.hdr["TELRA"] == ""


def test_TELDEC(hdr):
    assert hdr.hdr["TELDEC"] == ""


def test_RA(hdr):
    assert hdr.hdr["RA"] == ""


def test_DEC(hdr):
    assert hdr.hdr["DEC"] == ""


def test_RA_DEG(hdr):
    assert hdr.hdr["RA_DEG"] == ""


def test_DEC_DEG(hdr):
    assert hdr.hdr["DEC_DEG"] == ""


def test_EPOCH(hdr):
    assert hdr.hdr["EPOCH"] == ""


def test_HA(hdr):
    assert hdr.hdr["HA"] == ""


def test_HAD(hdr):
    assert hdr.hdr["HAD"] == ""


def test_AIRMASS(hdr):
    assert hdr.hdr["AIRMASS"] == ""


def test_AMSTART(hdr):
    assert hdr.hdr["AMSTART"] == ""


def test_AMEND(hdr):
    assert hdr.hdr["AMEND"] == ""


def test_XOFFSET(hdr):
    assert hdr.hdr["XOFFSET"] == ""


def test_YOFFSET(hdr):
    assert hdr.hdr["YOFFSET"] == ""


def test_PA(hdr):
    assert hdr.hdr["PA"] == ""


def test_CAMFOC(hdr):
    assert hdr.hdr["CAMFOC"] == ""


def test_TELFOC(hdr):
    assert hdr.hdr["TELFOC"] == ""


def test_TCSUT(hdr):
    assert hdr.hdr["TCSUT"] == ""


def test_TCSST(hdr):
    assert hdr.hdr["TCSST"] == ""


def test_MOONDIST(hdr):
    assert hdr.hdr["MOONDIST"] == ""


def test_MOONALT(hdr):
    assert hdr.hdr["MOONALT"] == ""


def test_MOONPHAS(hdr):
    assert hdr.hdr["MOONPHAS"] == ""


def test_ADC(hdr):
    assert hdr.hdr["ADC"] == ""


def test_ADCSPEED(hdr):
    assert hdr.hdr["ADCSPEED"] == ""


def test_CALMIRR(hdr):
    assert hdr.hdr["CALMIRR"] == ""


def test_GUIDING(hdr):
    assert hdr.hdr["GUIDING"] == ""


def test_GUIDTEXP(hdr):
    assert hdr.hdr["GUIDTEXP"] == ""


def test_GUIDFREQ(hdr):
    assert hdr.hdr["GUIDFREQ"] == ""


def test_GUIDOBJX(hdr):
    assert hdr.hdr["GUIDOBJX"] == ""


def test_GUIDOBJY(hdr):
    assert hdr.hdr["GUIDOBJY"] == ""


def test_AVGYCORR(hdr):
    assert hdr.hdr["AVGYCORR"] == ""


def test_AVGXCORR(hdr):
    assert hdr.hdr["AVGXCORR"] == ""


def test_GFOCUS(hdr):
    assert hdr.hdr["GFOCUS"] == ""


def test_THARLAMP(hdr):
    assert hdr.hdr["THARLAMP"] == ""


def test_HALLAMP(hdr):
    assert hdr.hdr["HALLAMP"] == ""


def test_AGITATOR(hdr):
    assert hdr.hdr["AGITATOR"] == ""


def test_THARMIRR(hdr):
    assert hdr.hdr["THARMIRR"] == ""


def test_DSTATUS(hdr):
    assert hdr.hdr["DSTATUS"] == ""


def test_DPOS(hdr):
    assert hdr.hdr["DPOS"] == ""


def test_DTEMP(hdr):
    assert hdr.hdr["DTEMP"] == ""


def test_DHUM(hdr):
    assert hdr.hdr["DHUM"] == ""


def test_DLAMP(hdr):
    assert hdr.hdr["DLAMP"] == ""


def test_DFLAT(hdr):
    assert hdr.hdr["DFLAT"] == ""


def test_TEMPEXT(hdr):
    assert hdr.hdr["TEMPEXT"] == ""


def test_PRESSURE(hdr):
    assert hdr.hdr["PRESSURE"] == ""


def test_HUMIDITY(hdr):
    assert hdr.hdr["HUMIDITY"] == ""


def test_WINDDIR(hdr):
    assert hdr.hdr["WINDDIR"] == ""


def test_WINSPEED(hdr):
    assert hdr.hdr["WINSPEED"] == ""


def test_SHTTTL(hdr):
    assert hdr.hdr["SHTTTL"] == ""


def test_COOLER(hdr):
    assert hdr.hdr["COOLER"] == ""


def test_TEMP(hdr):
    assert hdr.hdr["TEMP"] == dic["ccd_temp"]


def test_TEMPST(hdr):
    assert hdr.hdr["TEMPST"] == ""


def test_SERN(hdr):
    assert hdr.hdr["SERN"] == serial_number


def test_IM_DT0(hdr):
    assert hdr.hdr["IM_DT0"] == ""


def test_OBSTITLE(hdr):
    assert hdr.hdr["OBSTITLE"] == ""
