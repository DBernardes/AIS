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

    def __init__(self, air_mass: int | float, sky_condition: str):
        """Initialize the class.

        Parameters
        ----------
        air_mass: int | float
            The air mass in the light path.

        sky_condition: [photometric, regular, good]
            The condition of the sky at the observation night.
        """

        self.air_mass = air_mass
        self.sky_condition = sky_condition

    def _read_spreadsheet(self):
        spreasheet = pd.read_csv(self._SPECTRAL_RESPONSE_FILE, dtype=np.float64)
        self.atm_wavelength_interval = spreasheet["wavelength"]
        self.extinction_coefs = spreasheet[self.sky_condition] / 100

    def _calculate_atmosphere_transmitance(self, wavelength_interval):
        transmitance = [10 ** (-0.4 * k * self.air_mass) for k in self.extinction_coefs]
        spl = splrep(self.atm_wavelength_interval, transmitance)
        transmitance = splev(wavelength_interval, spl)
        return transmitance

    def apply_atmosphere_spectral_response(
        self,
        specific_photons_per_second,
        wavelength_interval,
    ):
        """Apply the atmosphere spectral response.

        This function applies the atmosphere spectral response on the
        calculated star specific flux

        Parameters
        ----------
        specific_photons_per_second : array like
            Photons per second per wavelength of the star

        wavelength_interval: array like
            Array with the wavelengths of the spectrum of the star.

        Returns
        -------
        specific_photons_per_second : array like
            Photons per second per wavelength of the star after the application atmosphere
            spectral response.
        """

        self._read_spreadsheet()
        transmitance = self._calculate_atmosphere_transmitance(wavelength_interval)
        specific_photons_per_second[0] = np.multiply(
            specific_photons_per_second[0], transmitance
        )

        return specific_photons_per_second
