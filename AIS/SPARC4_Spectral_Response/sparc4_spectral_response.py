"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""
import numpy as np
import pandas as pd


class Abstract_SPARC4_Spectral_Response:
    """Abstract class of the SPARC4 spectral response."""

    _CHANNEL_ID = 0
    _DIR_PATH = "./SPARC4_Spectral_Response/"

    def __init__(self):
        """Initialize the class."""

    def get_channel_ID(self):
        """Return the chanel ID.

        Returns
        -------
        _CHANNEL_ID : [1, 2, 3, 4]
            Channel ID
        """
        return self._CHANNEL_ID

    def write_spectrum(self, spectrum):
        """Write spectrum.

        This function writes the spectrum of the object in the class.

        Parameters
        ----------

        spectrum: array-like
            Spectrum of the object.
        """

        temp = np.asarray(spectrum)
        self.spectrum_length = len(spectrum)
        spectrum = np.zeros((4, self.spectrum_length))
        spectrum[0, :] = temp
        self.spectrum = spectrum

    def get_spectrum(self):
        """Get spectrum.

        This function returns the spectrum of the object.
        """

        return self.spectrum

    def calibration_wheel(self):
        """Calibration wheel spectrum response.

        Apllies the calibration wheel spectral response on the flux.
        """
        file = self._DIR_PATH + "calibration_wheel.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def retarder(self):
        """Retarder spectrum response.

        Apllies the retarder spectral response on the flux.
        """
        file = self._DIR_PATH + "retarder.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def analyzer(self):
        """Analyzer spectrum response.

        Apllies the analyzer spectral response on the flux.
        """
        file = self._DIR_PATH + "analyser.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def collimator(self):
        """Collimator spectrum response.

        Apllies the collimator spectral response on the flux.
        """
        file = self._DIR_PATH + "collimator.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def dichroic(self):
        """Dichroic spectrum response.

        Apllies the dichroic spectral response on the flux.
        """
        file = self._DIR_PATH + f"Channel {self._CHANNEL_ID}/" + "dichroic.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def camera(self):
        """Camera spectrum response.

        Apllies the camera spectral response on the flux.
        """
        file = self._DIR_PATH + f"Channel {self._CHANNEL_ID}/" + "camera.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def ccd(self):
        """CCD spectrum response.

        Apllies the CCD spectral response on the flux.
        """
        file = self._DIR_PATH + f"Channel {self._CHANNEL_ID}/" + "ccd.xlsx"
        stokes = self._read_spreadsheet(file)

        for i in range(self.spectrum_length):
            self.spectrum[:, i] = np.dot(stokes, self.spectrum[:, i])

    def _read_spreadsheet(self, file):
        ss = np.asarray(pd.read_excel(file))
        return ss


class Concrete_SPARC4_Spectral_Response_1(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 1."""

    _CHANNEL_ID = 1

    pass


class Concrete_SPARC4_Spectral_Response_2(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 2."""

    _CHANNEL_ID = 2

    pass


class Concrete_SPARC4_Spectral_Response_3(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 3."""

    _CHANNEL_ID = 3
    pass


class Concrete_SPARC4_Spectral_Response_4(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 4."""

    _CHANNEL_ID = 4
    pass
