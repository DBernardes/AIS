# -*- coding: utf-8 -*-
"""Header class tests.

This script tests the operation of the Header Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


import datetime
import os

import astropy.io.fits as fits
import pandas as pd
import pytest

from AIS.Header import Header

ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 1024,
}
ccd_temp = -70


@pytest.fixture
def hdr():
    hdr = Header(ccd_operation_mode, ccd_temp, 1)
    return hdr


file = os.path.join("AIS", "Header", "header.csv")
ss = pd.read_csv(file, sep=";")
ccd_gain = 3.37


# ----------------------- Initilization -----------------------------


def test_ccd_operation_mode(hdr):
    assert hdr.ccd_operation_mode == ccd_operation_mode


def test_ccd_temp(hdr):
    assert hdr.ccd_temp == ccd_temp


def test_channel(hdr):
    assert hdr.channel == 1


# ---------------------------- teste read spreadsheet ----------------------


def test_read_spreadsheet(hdr):
    hdr._read_spreadsheet(file)


def test_get_ccd_gain(hdr):
    assert hdr.ccd_gain == ccd_gain


# -----------------------------test_create_image_header---------------------

cards = [
    (keyword, "", comment) for keyword, comment in zip(ss["Keyword"], ss["Comment"])
]
dic = ccd_operation_mode
channel = 1
header = fits.Header(cards)
header["SIMPLE"] = True
header["NAXIS"] = 2

header["NAXIS1"] = dic["image_size"]
header["NAXIS2"] = dic["image_size"]
header["OBSERVER"] = "Johannes Kepler"
header["OBJECT"] = "HD5980"
header["INSTRUME"] = "SPARC4"
header["OBSTYPE"] = "NONE"
header["CCDSERN"] = channel + 9913
header["CHANNEL"] = channel
header["OBSLONG"] = -45.5825
header["OBSLAT"] = -22.53444444444445
header["OBSALT"] = 1864.0

date = datetime.datetime.now()
date_obs = date.strftime("%Y-%m-%dT%H:%M:%S")
header["DATE-OBS"] = date_obs
header["UTDATE"] = date_obs.split("T")[0]
header["UTTIME"] = date_obs.split("T")[1]

header["NCYCLES"] = 1
header["CYCLIND"] = 1
header["NFRAMES"] = 1
header["FRAMEIND"] = 1

header["EXPTIME"] = dic["t_exp"]
header["ACQMODE"] = "Kinetic"
header["PREAMP"] = "Gain " + str(dic["preamp"]) + "x"
header["READRATE"] = dic["readout"]
header["VSHIFT"] = 4.33
header["TRIGGER"] = "External"
header["EMMODE"] = dic["em_mode"]
header["EMGAIN"] = dic["em_gain"]
header["HBIN"] = dic["binn"]
header["VBIN"] = dic["binn"]
header["INITLIN"] = 1
header["INITCOL"] = 1
header["FINALLIN"] = dic["image_size"]
header["FINALCOL"] = dic["image_size"]
header["SHUTTER"] = "CLOSED"
header["COOLER"] = "ON"
header["CCDTEMP"] = ccd_temp
header["TGTEMP"] = ccd_temp
header["TEMPST"] = "TEMPERATURE_STABILIZED"
header["FRAMETRF"] = "ON"
header["VCLKAMP"] = "Normal"

header["GAIN"] = 3.37
header["RDNOISE"] = 6.66

header["RA"] = "00:00:00"
header["DEC"] = "00:00:00"
header["EQUINOX"] = 2000.0
header["TELFOCUS"] = 0.0
header["TCSHA"] = "00:00:00"
header["TCSDATE"] = date.strftime("%Y/%m/%d %H:%M:%S")

header["EXTTEMP"] = 11.7
header["AIRMASS"] = 1.01
header["PRESSURE"] = 760
header["HUMIDITY"] = 86.0

header["ACSVRSN"] = "1.0.0"
header["CTRLINTE"] = "S4GUI"
header["ACSMODE"] = "Real"
header["TCSMODE"] = "Real"


def test_create_header(hdr):
    new_header = hdr.create_header()
    del new_header["DATE-OBS"]
    del new_header["UTTIME"]
    del new_header["TCSDATE"]
    del header["DATE-OBS"]
    del header["UTTIME"]
    del header["TCSDATE"]
    assert header == new_header
