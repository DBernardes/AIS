"""
Header Class
============


This class creates the header that will be written in the images created by
the AIS.
"""

import os
from cmath import exp
import datetime

import astropy.io.fits as fits
import pandas as pd


class Header:

    """Image header class.

    Parameters
    ----------
    ccd_operation_mode : dictionary
        Operation mode parameters of the SPARC4 camera

        em_mode : ['EM', 'Conv']
            Electron Multiplying mode of the camera.
        em_gain : float
            EM gain of the camera.
        readout : [30, 20, 10, 1, 0.1]
            Readout rate of the pixels in MHz.
        preamp : [1, 2]
            Pre-amplifier gain
        binn : [1, 2]
            Binning of the pixels
        t_exp : float
            Exposure time in seconds        
        image_size : int
            Image size in pixels
    ccd_temp: float
        The CCD temperatura in celsius degree.    
    channel:
        The channel related to the camera.
    """

    _CSV_HEADER_FILE = os.path.join("AIS", "Header", "header.csv")
    _CSV_GAINS_FILE = os.path.join("AIS", "Header", "preamp_gains.csv")

    def __init__(self, ccd_operation_mode: dict, ccd_temp: float, channel: int):
        """Initialize the Header Class."""
        self.ccd_operation_mode = ccd_operation_mode
        self.ccd_temp = ccd_temp
        self.channel = channel

        self.get_ccd_gain()
        ss = self._read_spreadsheet(self._CSV_HEADER_FILE)
        cards = [(keyword, value, comment) for keyword, value,
                 comment in zip(ss['keyword'], ss['value'], ss['comment'])]
        self.header = fits.Header(cards)

    def get_ccd_gain(self):
        idx_tab = 0
        readout = self.ccd_operation_mode['readout']
        if self.ccd_operation_mode['em_mode'] == 'EM':
            idx_tab = [30, 20, 10, 1].index(readout) * 2
        else:
            idx_tab = [1, 0.1].index(readout) * 2 + 8
        idx_tab += self.ccd_operation_mode['preamp'] - 1
        ss = pd.read_csv(self._CSV_GAINS_FILE)
        self.ccd_gain = float(ss[f'{9913 + self.channel}'][idx_tab])

    @staticmethod
    def _read_spreadsheet(file):
        ss = pd.read_csv(file)
        return ss

    def create_header(self):
        """Create the image header.

        This functions returns a astropy.io.fits.Header class with all the information used 
        to create the artificial image.
        """
        dic = self.ccd_operation_mode
        self.header["ACQMODE"] = "Single"
        self.header["READMODE"] = "Image"
        self.header["HBIN"] = dic['binn']
        self.header["VBIN"] = dic['binn']
        self.header["TRIGGER"] = "External"
        self.header["EXPTIME"] = dic['t_exp']
        self.header["TEMP"] = self.ccd_temp
        self.header["READOUT"] = str(1 / dic['readout']) + "E-006"
        self.header["VSHIFT"] = "0.6E-06"
        self.header["GAIN"] = self.ccd_gain
        self.header["EMMODE"] = dic['em_mode']
        self.header["EMGAIN"] = dic['em_gain']
        self.header["PREAMP"] = str(dic['preamp']) + "x"
        self.header["SERN"] = self.channel + 9913
        self.header["CHANNEL"] = self.channel

        date = datetime.datetime.now()
        date_str = date.strftime('%Y%m%dT%H:%M:%S')
        self.header["DATE-OBS"] = date_str
        self.header["FILENAME"] = date_str

        return self.header
