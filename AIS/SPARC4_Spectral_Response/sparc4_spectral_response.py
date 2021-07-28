"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""


import os
import sys
from sys import exit

import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import splev, splrep


class Abstract_SPARC4_Spectral_Response:
    """Abstract class of the SPARC4 spectral response."""

    _CHANNEL_ID = 0
    _DIR_PATH = "SPARC4_Spectral_Response"

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

    def write_specific_flux(self, specific_flux, wavelength_interval):
        """Write the specific flux.

        This function writes the specific flux of the object in the class.

        Parameters
        ----------

        specific_flux: array-like
            Specific flux of the object.

        wavelength_interval: array-like.
            Wavelength interval of the specific flux.
        """
        self.specific_flux = specific_flux
        self.wavelength_interval = wavelength_interval
        self.specific_flux_length = len(specific_flux)

    def get_specific_flux(self):
        """Get the specific flux.

        Returns
        ------

        speficif_flux: array-like
            The specific flux of the object.
        """

        return self.specific_flux

    def get_specific_ordinary_ray(self):
        """Get the ordinary ray."""
        return self.specific_ordinary_ray

    def get_specific_extra_ordinary_ray(self):
        """Get the extra ordinary ray."""
        return self.specific_extra_ordinary_ray

    def apply_calibration_wheel(self):
        """Apply calibration wheel spectral response."""

        file = os.path.join(self._DIR_PATH, "calibration_wheel.xlsx")
        stokes = np.asarray(pd.read_excel(file))
        self.specific_flux = self._multiply_matrices(stokes, self.specific_flux)

    def apply_retarder(self):
        """Apply retarder spectral response."""

        file = os.path.join(self._DIR_PATH, "retarder.xlsx")
        stokes = np.asarray(pd.read_excel(file))
        self.specific_flux = self._multiply_matrices(stokes, self.specific_flux)

    def apply_analyser(self):
        """Apply analyser spectral response."""

        file = os.path.join(self._DIR_PATH, "analyser_ordinary.xlsx")
        stokes = np.asarray(pd.read_excel(file))
        self.specific_ordinary_ray = self._multiply_matrices(
            stokes, self.specific_flux.copy()
        )

        file = os.path.join(self._DIR_PATH, "analyser_extra_ordinary.xlsx")
        stokes = np.asarray(pd.read_excel(file))
        self.specific_extra_ordinary_ray = self._multiply_matrices(
            stokes, self.specific_flux.copy()
        )

    def apply_collimator(self):
        """Collimator spectral response."""

        file = os.path.join(self._DIR_PATH, "collimator.xlsx")
        coll_wavelength_interv, coll_transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(coll_transmitance, coll_wavelength_interv)

        try:
            self.specific_ordinary_ray = np.multiply(
                self.specific_ordinary_ray[0, :], transmitance
            )
            self.specific_extra_ordinary_ray = np.multiply(
                self.specific_extra_ordinary_ray[0, :], transmitance
            )
        except Exception:
            self.specific_ordinary_ray = np.multiply(
                self.specific_flux[0, :], transmitance
            )
            self.specific_extra_ordinary_ray = 0

    def apply_dichroic(self):
        """Apply the dichroic spectral response.


        This functions applies the spectral response of the two
        dichroics that compose each channel."""

        file = os.path.join(
            self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "dichroic 1.xlsx"
        )
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_ordinary_ray = np.multiply(
            self.specific_ordinary_ray, transmitance
        )
        self.specific_extra_ordinary_ray = np.multiply(
            self.specific_extra_ordinary_ray, transmitance
        )

        file = os.path.join(
            self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "dichroic 2.xlsx"
        )
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_ordinary_ray = np.multiply(
            self.specific_ordinary_ray, transmitance
        )
        self.specific_extra_ordinary_ray = np.multiply(
            self.specific_extra_ordinary_ray, transmitance
        )

    def apply_camera(self):
        """Apply the camera spectral response."""

        file = os.path.join(
            self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "camera.xlsx"
        )
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_ordinary_ray = np.multiply(
            self.specific_ordinary_ray, transmitance
        )
        self.specific_extra_ordinary_ray = np.multiply(
            self.specific_extra_ordinary_ray, transmitance
        )

    def apply_ccd(self):
        """Apply the ccd spectral response."""

        file = os.path.join(self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "ccd.xlsx")
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_ordinary_ray = np.multiply(
            self.specific_ordinary_ray, transmitance
        )
        self.specific_extra_ordinary_ray = np.multiply(
            self.specific_extra_ordinary_ray, transmitance
        )

    def _read_spreadsheet(self, file):
        ss = np.asarray(pd.read_excel(file))
        wavelength = [float(value) for value in ss[1:, 0]]
        transmitance = [float(value) / 100 for value in ss[1:, 1]]
        return wavelength, transmitance

    def _multiply_matrices(self, stokes, specific_flux):
        for i in range(len(specific_flux[0])):
            specific_flux[:, i] = np.dot(stokes, specific_flux[:, i])
        return specific_flux

    def _calculate_spline(self, transmitance, component_wavelength_interv):
        spl = splrep(component_wavelength_interv, transmitance)
        transmitance = splev(self.wavelength_interval, spl)
        return transmitance


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
