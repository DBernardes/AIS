import datetime
import os
from cmath import exp

import astropy.io.fits as fits
import pandas as pd

__all__ = ["Header"]


class Header:
    _local_path = os.path.dirname(__file__)
    _CSV_HEADER_FILE = os.path.join(_local_path, "header.csv")
    _CSV_GAINS_FILE = os.path.join(_local_path, "preamp_gains.csv")
    _CSV_READNOISE_FILE = os.path.join(_local_path, "read_noises.csv")

    def __init__(self, ccd_operation_mode: dict, ccd_temp: float, channel: int) -> None:
        """Initialize the class.

        Parameters
        ----------
        ccd_operation_mode : dict
            Operation mode parameters of the SPARC4 camera.
        ccd_temp : float
            The CCD temperatura in celsius degree.
        channel : [1, 2, 3, or 4].
            The channel ID of the camera
        """
        self.ccd_operation_mode = ccd_operation_mode
        self.ccd_temp = ccd_temp
        self.channel = channel

        idx_tab = self._find_index_tab()
        ss = pd.read_csv(self._CSV_GAINS_FILE)
        self.ccd_gain = float(ss[f"{9913 + self.channel}"][idx_tab])
        ss = pd.read_csv(self._CSV_READNOISE_FILE)
        self.read_noise = float(ss[f"{9913 + self.channel}"][idx_tab])

        ss = self._read_spreadsheet(self._CSV_HEADER_FILE)
        cards = [
            (keyword, "", comment)
            for keyword, comment in zip(ss["Keyword"], ss["Comment"])
        ]
        self.header = fits.Header(cards)

        return

    def _find_index_tab(self) -> int:
        idx_tab = 0
        readout = self.ccd_operation_mode["readout"]
        if self.ccd_operation_mode["em_mode"] == "EM":
            idx_tab = [30, 20, 10, 1].index(readout) * 2
        else:
            idx_tab = [1, 0.1].index(readout) * 2 + 8
        idx_tab += self.ccd_operation_mode["preamp"] - 1
        return idx_tab

    @staticmethod
    def _read_spreadsheet(file) -> pd.DataFrame:
        ss = pd.read_csv(file, sep=";")
        return ss

    def create_header(self) -> fits.Header:
        """Create the image header.

        This functions returns a astropy.io.fits.Header class with all
        the information used to create the artificial image.
        """
        dic = self.ccd_operation_mode
        self.header["SIMPLE"] = True
        self.header["NAXIS"] = 2

        self.header["NAXIS1"] = dic["image_size"]
        self.header["NAXIS2"] = dic["image_size"]
        self.header["OBSERVER"] = "Johannes Kepler"
        self.header["OBJECT"] = "HD5980"
        self.header["INSTRUME"] = "SPARC4"
        self.header["OBSTYPE"] = "NONE"
        self.header["CCDSERN"] = self.channel + 9913
        self.header["CHANNEL"] = self.channel
        self.header["OBSLONG"] = -45.5825
        self.header["OBSLAT"] = -22.53444444444445
        self.header["OBSALT"] = 1864.0

        date = datetime.datetime.now()
        date_obs = date.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
        self.header["DATE-OBS"] = date_obs
        self.header["UTDATE"] = date_obs.split("T")[0]
        self.header["UTTIME"] = date_obs.split("T")[1]

        self.header["NCYCLES"] = 1
        self.header["CYCLIND"] = 1
        self.header["NFRAMES"] = 1
        self.header["FRAMEIND"] = 1

        self.header["EXPTIME"] = dic["t_exp"]
        self.header["ACQMODE"] = "Kinetic"
        self.header["PREAMP"] = "Gain " + str(dic["preamp"]) + "x"
        self.header["READRATE"] = dic["readout"]
        self.header["VSHIFT"] = 4.33
        self.header["TRIGGER"] = "External"
        self.header["EMMODE"] = dic["em_mode"]
        self.header["EMGAIN"] = dic["em_gain"]
        self.header["HBIN"] = dic["binn"]
        self.header["VBIN"] = dic["binn"]
        self.header["INITLIN"] = 1
        self.header["INITCOL"] = 1
        self.header["FINALLIN"] = dic["image_size"]
        self.header["FINALCOL"] = dic["image_size"]
        self.header["SHUTTER"] = "CLOSED"
        self.header["COOLER"] = "ON"
        self.header["CCDTEMP"] = self.ccd_temp
        self.header["TGTEMP"] = self.ccd_temp
        self.header["TEMPST"] = "TEMPERATURE_STABILIZED"
        self.header["FRAMETRF"] = "ON"
        self.header["VCLKAMP"] = "Normal"

        self.header["GAIN"] = self.ccd_gain
        self.header["RDNOISE"] = self.read_noise

        self.header["RA"] = "00:00:00"
        self.header["DEC"] = "00:00:00"
        self.header["EQUINOX"] = 2000.0
        self.header["TELFOCUS"] = 0.0
        self.header["TCSHA"] = "00:00:00"
        self.header["TCSDATE"] = date.strftime("%Y/%m/%d %H:%M:%S")

        self.header["EXTTEMP"] = 11.7
        self.header["AIRMASS"] = 1.01
        self.header["PRESSURE"] = 760
        self.header["HUMIDITY"] = 86.0

        self.header["ACSVRSN"] = "1.0.0"
        self.header["CTRLINTE"] = "S4GUI"
        self.header["ACSMODE"] = "Real"
        self.header["TCSMODE"] = "Real"

        return self.header
