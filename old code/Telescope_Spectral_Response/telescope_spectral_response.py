"""
Telescope Spectral Response Class
===================================

This class calculates the output flux of an astronomical object as a funtion
of the 1.6 m Perkin-Elmer spectral response.
"""

import os

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep


class Telescope_Spectral_Response:

    """Telescope Spectral Response Class."""

    _SPECTRAL_RESPONSE_FILE = os.path.join(
        "AIS", "Telescope_Spectral_Response", "telescope_spectral_response.csv"
    )

    def __init__(self):
        """Initilize the class."""
        pass

    def _read_spreadsheet(self):
        ss = pd.read_csv(
            self._SPECTRAL_RESPONSE_FILE, dtype=np.float64, skiprows=1, decimal="."
        )
        self.tel_wavelength_interval = ss["(nm)"]
        self.reflectance = ss["(%)"] / 100

    def _calculate_spline(self, wavelength_interval):
        spl = splrep(self.tel_wavelength_interval, self.reflectance)
        reflectance = splev(wavelength_interval, spl)
        return reflectance

    def apply_telescope_spectral_response(
        self, star_specific_photons_per_second, wavelength_interval
    ):
        """
        Apply the telescope spectral response.

        This function applies the telescope spectral response on the
        calculated star specific flux.

        Parameters
        ----------
        star_specific_photons_per_second : array like
            Specific flux of the star

        wavelength_interval: array like
            Array with the wavelengths of the spectrum of the star.

        Returns
        -------
        star_specific_photons_per_second : array like
            Specific flux of the star generated by the spectral response of the
            1.6 m Perkin-Elmer telescope
        """
        self.star_specific_photons_per_second = star_specific_photons_per_second
        self._read_spreadsheet()
        reflectance = self._calculate_spline(wavelength_interval)
        new_star_specific_photons_per_second = np.multiply(
            star_specific_photons_per_second[0], reflectance
        )
        self.star_specific_photons_per_second[0] = new_star_specific_photons_per_second

        return star_specific_photons_per_second