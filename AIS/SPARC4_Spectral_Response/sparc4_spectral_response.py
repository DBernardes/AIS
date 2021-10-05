"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""


import os

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep

# from sys import exit


class Abstract_SPARC4_Spectral_Response:

    """Abstract class of the SPARC4 spectral response."""

    _CHANNEL_ID = 0
    _DIR_PATH = "SPARC4_Spectral_Response"

    def __init__(self, wavelength_interval):
        """
        Initialize the class.

        Parameters
        ----------

        wavelength_interval: array-like.
            Wavelength interval of the specific flux.
        """

        self.wavelength_interval_len = len(wavelength_interval)
        self.wavelength_interval = wavelength_interval

    def get_channel_ID(self):
        """
        Return the chanel ID.

        Returns
        -------
        _CHANNEL_ID : [1, 2, 3, 4]
            Channel ID
        """
        return self._CHANNEL_ID

    def write_specific_flux(self, specific_flux):
        """
        Write the specific flux.

        This function writes the specific flux of the object in the class.

        Parameters
        ----------

        specific_flux: array-like
            Specific flux of the object.
        """
        self.specific_ordinary_ray = specific_flux
        self.specific_extra_ordinary_ray = np.zeros((4, self.wavelength_interval_len))

    def get_specific_ordinary_ray(self):
        """Get the ordinary ray."""
        return self.specific_ordinary_ray

    def get_specific_extra_ordinary_ray(self):
        """Get the extra ordinary ray."""
        return self.specific_extra_ordinary_ray

    def apply_calibration_wheel(self):
        """Apply calibration wheel spectral response."""
        file = os.path.join(self._DIR_PATH, "calibration_wheel.csv")
        matrix = np.loadtxt(open(file, "rb"), delimiter=",")
        self.specific_ordinary_ray = self._multiply_matrices(
            matrix, self.specific_ordinary_ray
        )

    def apply_retarder(self):
        """Apply retarder spectral response."""
        file = os.path.join(self._DIR_PATH, "retarder.csv")
        matrix = np.loadtxt(open(file, "rb"), delimiter=",")
        self.specific_ordinary_ray = self._multiply_matrices(
            matrix, self.specific_ordinary_ray
        )

    def apply_analyser(self):
        """Apply analyser spectral response."""
        file = os.path.join(self._DIR_PATH, "analyser_ordinary.csv")
        matrix = np.loadtxt(open(file, "rb"), delimiter=",")
        temp = self.specific_ordinary_ray
        self.specific_ordinary_ray = self._multiply_matrices(matrix, temp.copy())

        file = os.path.join(self._DIR_PATH, "analyser_extra_ordinary.csv")
        matrix = np.loadtxt(open(file, "rb"), delimiter=",")
        self.specific_extra_ordinary_ray = self._multiply_matrices(matrix, temp.copy())

    def apply_collimator(self):
        """Collimator spectral response."""
        file = os.path.join(self._DIR_PATH, "collimator.csv")
        self._apply_optical_component(file)

    def apply_dichroic(self):
        """
        Apply the dichroic spectral response.


        This functions applies the spectral response of the two
        dichroics that compose each channel.
        """
        for i in [1, 2]:
            file = os.path.join(
                self._DIR_PATH, f"Channel {self._CHANNEL_ID}", f"dichroic_{i}.csv"
            )
            self._apply_optical_component(file)

    def apply_camera(self):
        """Apply the camera spectral response."""
        file = os.path.join(self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "camera.csv")
        self._apply_optical_component(file)

    def apply_ccd(self):
        """Apply the ccd spectral response."""
        file = os.path.join(self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "ccd.csv")
        self._apply_optical_component(file)

    def _read_spreadsheet(self, file):
        ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
        wavelength = ss["(nm)"]
        transmitance = ss["(%)"] / 100
        return wavelength, transmitance

    def _multiply_matrices(self, matrix, specific_flux):
        for i in range(len(specific_flux[0])):
            specific_flux[:, i] = np.dot(matrix, specific_flux[:, i])
        return specific_flux

    def _calculate_spline(self, transmitance, component_wavelength_interv):
        spl = splrep(component_wavelength_interv, transmitance)
        transmitance = splev(self.wavelength_interval, spl)
        return transmitance

    def _apply_optical_component(self, file):
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_ordinary_ray = np.multiply(
            self.specific_ordinary_ray, transmitance
        )
        self.specific_extra_ordinary_ray = np.multiply(
            self.specific_extra_ordinary_ray, transmitance
        )


class Concrete_SPARC4_Spectral_Response_1(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 1."""

    _CHANNEL_ID = 1


class Concrete_SPARC4_Spectral_Response_2(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 2."""

    _CHANNEL_ID = 2


class Concrete_SPARC4_Spectral_Response_3(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 3."""

    _CHANNEL_ID = 3


class Concrete_SPARC4_Spectral_Response_4(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 4."""

    _CHANNEL_ID = 4
