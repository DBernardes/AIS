"""
Header Class
============


This class creates the header that will be written in the images created by
the AIS.
"""

import os
from sys import exit

import astropy.io.fits as fits
import pandas as pd


class Header:
    """Image header class.

    Parameters
    ----------
    dic : dictionary
        Operation mode parameters of the SPARC4 camera

        em_mode : [0, 1]
            CCD Electron Multiplying Mode
        em_gain : float
            CCD Electron Multiplying gain
        hss : [0.1, 1, 10, 20, 30]
            Horizontal Shift Speed of the pixels
        preamp : [1, 2]
            Pre-amplifier gain
        binn : [1, 2]
            Binning of the pixels
        t_exp : float
            Exposure time in seconds
        ccd_temp : float
            CCD temperature in ºC
        image_size : int
            Image size in pixels
    ccd_gain : float
        CCD gain in e-/ADU
    serial_number:
        CCD serial number
    """

    _CSV_HEADER_FILE = os.path.join("Header", "header.csv")

    def __init__(self, ccd_operation_mode, ccd_gain, serial_number):
        """Initialize the Header Class."""
        self.em_mode = ccd_operation_mode["em_mode"]
        self.NOISE_FACTOR = 1
        self.em_gain = 1
        if self.em_mode == 1:
            self.NOISE_FACTOR = 1.41
            self.em_gain = ccd_operation_mode["em_gain"]
        self.preamp = ccd_operation_mode["preamp"]
        self.hss = ccd_operation_mode["hss"]
        self.binn = ccd_operation_mode["binn"]
        self.t_exp = ccd_operation_mode["t_exp"]
        self.ccd_temp = ccd_operation_mode["ccd_temp"]
        self.image_size = ccd_operation_mode["image_size"]
        self.ccd_gain = ccd_gain
        self.serial_number = serial_number

    def _read_spreadsheet(self):
        ss = pd.read_csv(self._CSV_HEADER_FILE)
        self.keywords = ss["Keyword"]
        self.comments = ss["Comment"]

    def create_header(self):
        """Create the image header.

        This functions writes a simple header with the used parameters for
        the CCD operation mode for the image FITS file
        """
        self._read_spreadsheet()
        n = len(self.keywords)
        self.hdr = fits.Header()
        self.hdr["SIMPLE"] = "T"
        self.hdr["BITPIX"] = "16"
        self.hdr["NAXIS1"] = self.image_size
        self.hdr["NAXIS2"] = self.image_size
        for i in range(n):
            self.hdr[self.keywords[i]] = ("", self.comments[i])
        self._write_header_values()
        return self.hdr

    def _write_header_values(self):
        """Write header values.

        This functions writes the values used to create the
        image in its respecive keyword.
        """
        self.hdr["EXTEND"] = "T"
        self.hdr["COMMENT"] = (
            "and Astrophysics, volume 376, page 359; bibcode:" + "2001A&A...376..3"
        )
        self.hdr["ACQMODE"] = "Single"
        self.hdr["READMODE"] = "Image"
        self.hdr["IMGRECT"] = f"1, {self.image_size}, {self.image_size}, 1"
        self.hdr["HBIN"] = self.binn
        self.hdr["VBIN"] = self.binn
        self.hdr["TRIGGER"] = "External"
        self.hdr["EXPTIME"] = self.t_exp
        self.hdr["TEMP"] = self.ccd_temp
        self.hdr["READOUT"] = str(1 / self.hss) + "E-006"
        self.hdr["VSHIFT"] = "4.33E-06"
        self.hdr["GAIN"] = self.ccd_gain
        em_mode = "Conventional"
        if self.em_mode == 1:
            em_mode = "Electron Multiplying"
        self.hdr["EMMODE"] = em_mode
        self.hdr["EMGAIN"] = self.em_gain
        self.hdr["PREAMP"] = str(self.preamp) + "x"
        self.hdr["SERN"] = self.serial_number
        self.hdr["DATE"] = "2017-07-14T00:00:58"
        self.hdr["FILENAME"] = "hats-24_I_transito_001"
