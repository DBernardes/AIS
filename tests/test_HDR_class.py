# -*- coding: utf-8 -*-
"""Header class tests.

This script tests the operation of the Header Class.

Created on Fri Apr 16 11:53:12 2021

@author: denis
"""


from ast import keyword
import astropy.io.fits as fits
import os
import datetime
import numpy as np
import pandas as pd
import pytest
from AIS.Header import Header
from tests.AIS_spectral_response_curves import ccd_operation_mode, ccd_temp


@pytest.fixture
def hdr():
    hdr = Header(ccd_operation_mode, ccd_temp, 1)
    return hdr


file = os.path.join("AIS", "Header", "header.csv")
ss = pd.read_csv(file, sep='\t')
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

cards = [(keyword, '', comment)
         for keyword, comment in zip(ss['keyword'], ss['comment'])]
dic = ccd_operation_mode
channel = 1
header = fits.Header(cards)
header["ACQMODE"] = "Single"
header["READMODE"] = "Image"
header["HBIN"] = dic['binn']
header["VBIN"] = dic['binn']
header["TRIGGER"] = "External"
header["EXPTIME"] = dic['t_exp']
header["TEMP"] = ccd_temp
header["READOUT"] = dic['readout']
header["VSHIFT"] = 4.33
header["GAIN"] = ccd_gain
header["EMMODE"] = dic['em_mode']
header["EMGAIN"] = dic['em_gain']
header["PREAMP"] = str(dic['preamp']) + "x"
header["SERN"] = 9913 + channel
header['CHANNEL'] = channel
date = datetime.datetime.now()

header["DATE-OBS"] = date.strftime('%Y%m%dT%H:%M:%S')
header["FILENAME"] = date.strftime('%Y%m%d') + f'_s4c{channel}_000000.fits'
header["UTDATE"] = date.strftime('%Y%m%d')
header["UTTIME"] = date.strftime('%H:%M:%S')


# def test_create_header(hdr):
#     new_header = hdr.create_header()
#     assert header == new_header
