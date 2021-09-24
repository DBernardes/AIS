"""
Atmosphere Spectral Response
============================

This class calculates the output flux of an astronomical object as a funtction
of the atmosphere spectral response.
"""
import os

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep


class Atmosphere_Spectral_Response:
    """Atmosphere Spectral Response Class."""

    _SPECTRAL_RESPONSE_FILE = os.path.join(
        "Atmosphere_Spectral_Response", "atmosphere_spectral_response.csv"
    )

    def __init__(self):
        """Initialize the class."""
        pass

    def _read_spreadsheet(self):
        ss = pd.read_csv(
            self._SPECTRAL_RESPONSE_FILE, dtype=np.float64, skiprows=1, decimal=","
        )
        self.atm_wavelength_interval = ss["(nm)"]
        self.transmitance = ss["(%)"] / 100

    def _calculate_spline(self, wavelength_interval):
        spl = splrep(self.atm_wavelength_interval, self.transmitance)
        transmitance = splev(wavelength_interval, spl)
        return transmitance

    def apply_atmosphere_spectral_response(
        self, star_specific_flux, l_init, l_final, l_step
    ):
        """Apply the atmosphere spectral response.

        This function applies the atmosphere spectral response on the
        calculated star specific flux

        Parameters
        ----------
        star_specific_flux : array like
            Specific flux of the star

        Returns
        -------
        star_specific_flux : array like
            Specific flux of the star after the application atmosphere
            spectral response.

        l_init: int
            Initial wavelength in nanometers.

        l_final: int
            Final wavelength in nanometers.

        l_step: int
            Step for the wavelength interval in nanometers.
        """
        self.star_specific_flux = star_specific_flux
        wavelength_interval = range(l_init, l_final, l_step)
        self._read_spreadsheet()
        transmitance = self._calculate_spline(wavelength_interval)
        new_specific_flux = np.multiply(star_specific_flux[0, :], transmitance)
        self.star_specific_flux[0, :] = new_specific_flux

        return self.star_specific_flux
