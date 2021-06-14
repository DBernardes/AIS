"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""


import sys

import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import splev, splrep


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

        This function returns the specific flux of the object.
        """

        return self.specific_flux

    def apply_photometric_component_spectral_response(self, name):
        """Apply photometric spectral response.

        Apllies the spectral response of a photometric component.

        Parameters
        ----------

        name: string
            The name of the photometric component.
        """
        file = self._DIR_PATH + f"Channel {self._CHANNEL_ID}/" + name + ".xlsx"
        wavelength_interv, transmitance = self._read_spreadsheet(file)
        new_transmitance = self._calculate_spline(transmitance, wavelength_interv)
        self.specific_flux = np.multiply(self.specific_flux, new_transmitance)

    def apply_polarimetric_component_spectral_response(self, name):
        """Apply polarimetric spectral response.

        Apllies the spectral response of the polarimetric component.

        Parameters
        ----------

        name: string
            The name of the polarimetric compoent.
        """

        file = self._DIR_PATH + name + ".xlsx"
        stokes = np.asarray(pd.read_excel(file))
        self._multiply_matrices(stokes, self.specific_flux)

    def collimator(self):
        """Collimator spectral response.

        Apllies the collimator spectral response on the flux.
        """
        file = self._DIR_PATH + "collimator.xlsx"
        coll_wavelength_interv, coll_transmitance = self._read_spreadsheet(file)
        transmitance = self._calculate_spline(coll_transmitance, coll_wavelength_interv)
        self.specific_flux = np.multiply(self.specific_flux[0, :], transmitance)

    def _read_spreadsheet(self, file):
        ss = np.asarray(pd.read_excel(file))
        wavelength = [float(value) for value in ss[1:, 0]]
        transmitance = [float(value) / 100 for value in ss[1:, 1]]
        return wavelength, transmitance

    def _multiply_matrices(self, a, b):
        for i in range(self.specific_flux_length):
            self.specific_flux[:, i] = np.dot(a, b[:, i])

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
