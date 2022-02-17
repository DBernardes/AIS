"""

Atmosphere Spectral Response
============================

This class calculates the output flux of an astronomical object as a funtction
of the atmosphere spectral response.
"""
import os
from sys import exit

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep


class Atmosphere_Spectral_Response:

    """Atmosphere Spectral Response Class."""

    _SPECTRAL_RESPONSE_FILE = os.path.join(
        "AIS", "Atmosphere_Spectral_Response", "atmosphere_spectral_response.csv"
    )

    def __init__(self, air_mass, sky_condition):
        """Initialize the class.

        Parameters
        ----------
        air_mass: float
            The air mass in the light path.

        sky_condition: {photometric, regular, good}
            The condition of the sky at the observation night.
        """

        self.air_mass = air_mass
        self.sky_condition = sky_condition

    def _read_spreadsheet(self):
        spreasheet = pd.read_csv(self._SPECTRAL_RESPONSE_FILE)
        self.atm_wavelength_interval = [
            float(value) for value in spreasheet["Wavelength"][1:]
        ]
        self.extinction_coefs = [
            float(value) / 100 for value in spreasheet[self.sky_condition][1:]
        ]

    def _calculate_atmosphere_transmitance(self, wavelength_interval):
        transmitance = [10 ** (-0.4 * k * self.air_mass) for k in self.extinction_coefs]
        spl = splrep(self.atm_wavelength_interval, transmitance)
        transmitance = splev(wavelength_interval, spl)
        return transmitance

    def apply_atmosphere_spectral_response(
        self,
        star_specific_flux,
        wavelength_interval,
    ):
        """Apply the atmosphere spectral response.

        This function applies the atmosphere spectral response on the
        calculated star specific flux

        Parameters
        ----------
        star_specific_flux : array like
            Specific flux of the star

        wavelength_interval: array like
            Array with the wavelengths of the spectrum of the star.

        Returns
        -------
        star_specific_flux : array like
            Specific flux of the star after the application atmosphere
            spectral response.
        """
        self.star_specific_flux = star_specific_flux
        self._read_spreadsheet()
        transmitance = self._calculate_atmosphere_transmitance(wavelength_interval)
        new_specific_flux = np.multiply(star_specific_flux[0, :], transmitance)
        self.star_specific_flux[0, :] = new_specific_flux

        return self.star_specific_flux
