"""
Spectrum Calculation Class
===========================

This class calculates the spectrum of an astronomical object as a function of its
magnitude
"""

import os
from importlib.util import spec_from_file_location
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep


class Spectrum_Calculation:
    """
    Spectrum Calculation class.

    This class calculates the star and the sky spectrum based on the object
    magnitude.
    """

    PLANCK_CONSTANT = 6.62607004e-34  # m2 kg / s
    LIGHT_SPPED = 3e8  # m/s
    TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    BAND_PASS = 0.2e-6  # m
    S_0 = 4e-2  # W/m2/m
    SPREADSHEET_MOON_MAGNITUDE = os.path.join(
        "AIS", "Spectrum_Calculation", "moon_magnitude.csv"
    )

    def __init__(
        self,
        wavelength_interval: Iterable,
    ):
        """
        Initialize the class.

        Parameters
        ----------

        wavelength_interval: array like
            Array with the wavelengths of the spectrum of the star.

        """
        self.wavelength_interval = wavelength_interval

    def _read_spreadsheet(self, moon_condition):
        spreasheet = pd.read_csv(self.SPREADSHEET_MOON_MAGNITUDE, dtype=np.float64)
        moon_magnitudes = spreasheet[moon_condition]
        moon_wavelenght = spreasheet["wavelength"]

        return moon_wavelenght, moon_magnitudes

    def _calculate_spline(self, component_wavelength_interv, value):
        spl = splrep(component_wavelength_interv, value)
        interpolated_value = splev(self.wavelength_interval, spl)
        return interpolated_value

    def calculate_sky_specific_photons_per_second(self, moon_condition: str):
        """Calculate the specific photons per second of the sky.

        Returns
        -------

        sky_specific_photons_per_second: array-like
            Specific phontons per second of the sky.
        """
        moon_wavelength, moon_magnitudes = self._read_spreadsheet(moon_condition)
        moon_magnitudes = self._calculate_spline(moon_wavelength, moon_magnitudes)

        temp = []
        for magnitude, wavelength in zip(moon_magnitudes, self.wavelength_interval):
            specific_photons_per_second = self._convert_magnitude(wavelength, magnitude)
            temp.append(specific_photons_per_second)

        sky_specific_photons_per_second = np.zeros((4, len(self.wavelength_interval)))
        sky_specific_photons_per_second[0, :] = temp

        return sky_specific_photons_per_second

    def _convert_magnitude(self, wavelength, magnitude):
        """Convert magnitude in photons per second, as a function of the wavelength.

        The calculation of the specific photons per second :math:`P_{\lambda}`
        is based on the expression

        .. math::
            P_{\lambda} = \frac{S_0 \times 10^{-mag/2.5} * \lambda * B * tel_{area}}{h * c}

        where :math:`S_0 = 4 \tim 10^{-2}` W/m2/m, :math:`mag` is the magnitude of the object,
        :math:`\lambda` is the wavelength in nm, :math:`B_P` is the width of the band pass being
        considered, :math:`tel_{area}` is the telescope effective area, :math:`h` is the Planck
        constant, and :math:`c` is the light speed.
        """
        wavelength *= 1e-9
        photons_number = (
            self.S_0
            * 10 ** (-magnitude / 2.5)
            * wavelength
            * self.BAND_PASS
            * self.TELESCOPE_EFFECTIVE_AREA
            / (self.PLANCK_CONSTANT * self.LIGHT_SPPED)
        )
        return photons_number

    def calculate_star_specific_photons_per_second(self, magnitude: int | float):
        """Calculate the specific photons per second of the star.

        Returns
        -------

        star_specific_photons_per_second: array-like
            Specific phontons per second of the star.
        """

        temp = []
        for wavelength in self.wavelength_interval:
            specific_photons_per_second = self._convert_magnitude(wavelength, magnitude)
            temp.append(specific_photons_per_second)

        star_specific_photons_per_second = np.zeros((4, len(self.wavelength_interval)))
        star_specific_photons_per_second[0, :] = temp

        return star_specific_photons_per_second
