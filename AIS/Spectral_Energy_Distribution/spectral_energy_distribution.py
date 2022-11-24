"""
Spectral Energy Distribution Class
==================================


The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
"""

from numpy import ndarray
import numpy as np
from scipy.interpolate import splev, splrep
import pandas as pd
from scipy.constants import c, h, k
from sbpy.calib import vega_fluxd
from math import pi
import os
from sys import exit

__all__ = ['Source', 'Sky']


class Spectral_Energy_Distribution:
    """Spectral Energy Distribution Class

    The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
    """
    EFFECT_WAVELENGTH = 555.6  # nm
    TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    S_0 = vega_fluxd.get()["Johnson V"].value*1e7  # W/(m.m2)
    BASE_PATH = os.path.join(
        'AIS', 'Spectral_Energy_Distribution')

    def __init__(self):
        """Initialize the Spectral Energy Distribution Class."""
        self.sed = np.linspace(100, 1000, 100)
        return

    def get_sed(self) -> ndarray:
        """Get the Spectral Energy Distribution."""
        return self.sed

    def write_sed(self, sed: ndarray):
        """Write the Spectral Energy Distribution to the class."""
        self.sed = sed

    def calculate_sed(self):
        """Calculate the Spectral Energy Distribution.

        This method is an abstract method that must be implemented in the child classes.

        Returns
        -------
        ndarray
            The Spectral Energy Distribution of the object in photons/s/m.
        """
        return np.linspace(100, 1000, 100)

    @staticmethod
    def _interpolate_spectral_distribution(wavelength, spectral_distribution, obj_wavelength):
        spl = splrep(wavelength, spectral_distribution)
        interpolated_spectral_distribution = splev(obj_wavelength, spl)
        return interpolated_spectral_distribution

    def _calculate_photons_density(self, magnitude) -> float:
        return self.S_0*10**(-magnitude/2.5)*self.TELESCOPE_EFFECTIVE_AREA*self.EFFECT_WAVELENGTH*1e-9/(h*c)


class Source(Spectral_Energy_Distribution):
    """Source Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the astronomical object.
    """

    def __init__(self):
        self.SPECTRAL_LIB_PATH = os.path.join(
            self.BASE_PATH, 'Spectral_Library')
        return

    def calculate_sed(self,
                      calculation_method: str,
                      magnitude: int | float,
                      wavelength_interval: tuple = (),
                      temperature: int | float = 0,
                      spectral_type: str = '') -> ndarray:
        """Get the Spectral Energy Distribution of the astronomical object.

        Parameters
        ----------

        calculation_method : ['blackbody', 'spectral_library']
            The method used to calculate the SED.
            If the user wants to use a blackbody SED, the calculation_method must be 'blackbody'.
            If the user wants to use a spectral standard SED, the calculation_method must be 'spectral_library'.

            In the 'blackbody' case, the spectral response of the object is calculated using the Planck function,
            given the temperature and the wavelength interval of the object. In the 'spectral_library' case, the
            spectral response and the wavelength of the object are obtained using a library of spectral types. 
            These spectrums are taken from the Library of Stellar Spectrum of the ESO, and they can be found at:
            https://www.eso.org/sci/observing/tools/standards/spectra/index.html.
            The level of the spectral response is adjusted using the magnitude of the object in the V band.

        magnitude : int | float
            The magnitude of the astronomical object in the V band.
            The magnitude is used to calculate the effective flux of the astronomical object.

        wavelength_interval : tuple, optional
            The wavelength interval, in nm, of the astronomical object.
            This parameter must be a tuple with three elements, where the first element is the initial wavelength,
            the second element is the final wavelength and the third element is the number of elements in the array.
            This parameter is used only if the calculation_method is 'blackbody'.

        temperature : int | float, optional
            The blackbody temperature of the astronomical object in Kelvin.
            This parameter is used only if the calculation_method is 'blackbody'.

        spectral_type : str, optional
            The spectral type of the star that will be used to calculate the SED.
            This parameter is used only if the calculation_method is 'spectral_standard'.
            The available spectral types can be found using the print_available_spectral_types() method.

        Returns
        -------
            wavelength : ndarray
                The wavelength of the astronomical object in nm.

            sed : ndarray
                The SED of the astronomical object in photons/m/s.
        """

        if calculation_method == 'blackbody':
            wv = wavelength_interval
            wavelength = np.linspace(wv[0], wv[1], wv[2])
            sed = self._calculate_sed_blackbody(
                wavelength, temperature)
            normalization_flux = self._interpolate_spectral_distribution(
                wavelength, sed, self.EFFECT_WAVELENGTH)
            sed /= normalization_flux
        elif calculation_method == 'spectral_library':
            wavelength, sed = self._read_spectral_library(
                spectral_type.lower())
        else:
            raise ValueError(
                "The calculation_method must be 'blackbody' or 'spectral_library'.")

        effective_flux = self._calculate_photons_density(
            magnitude)
        sed = sed * effective_flux
        return wavelength, sed

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature):
        return 8 * pi * h * c / (wavelength * 1e-9) ** 5 * 1 / (np.exp(h * c / (wavelength * 1e-9 * k * temperature)) - 1)

    def _read_spectral_library(self, spectral_type):
        path = os.path.join(self.SPECTRAL_LIB_PATH,
                            'uk' + spectral_type + '.csv')
        try:
            sed = pd.read_csv(path)
            return sed['wavelength (nm)'], sed['flux (F_lambda)']
        except FileNotFoundError:
            print(f"\nThe spectral type {spectral_type} is not available.")
            self.print_available_spectral_types()
            raise FileNotFoundError

    def print_available_spectral_types(self):
        """Print the available spectral types."""
        spec_types = os.listdir(self.SPECTRAL_LIB_PATH)
        print('\nAvailable spectral types:')
        print('-------------------------\n')
        spec_types = [spec_type.split('.')[0][2:] for spec_type in spec_types]
        print(*spec_types, sep=', ')


class Sky(Spectral_Energy_Distribution):
    """Sky Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the emission of the sky. 
    """
    CSV_FILE = 'moon_magnitude.csv'

    def __init__(self):
        '''Initialize the Sky class.'''
        return

    def _read_csv(self, file_name, value_name):
        file_name = os.path.join(self.BASE_PATH, file_name)
        ss = pd.read_csv(file_name)
        wavelenght = ss["wavelength"]
        value = ss[value_name]

        return wavelenght, value

    def calculate_sed(self, moon_phase: str, object_wavelength: ndarray) -> ndarray:
        """Get the Spectral Energy Distribution of the sky.

        Parameters
        ----------
        moon_phase : ['new', 'first quarter', 'third quarter', 'full']
            The phase of the moon.

        object_wavelength : ndarray
            The wavelength interval, in nm, of the astronomical object.

        Returns
        -------
            wavelength : ndarray
                The wavelength of the sky in nm.

            sed : ndarray
                The SED of the sky in photons/m/s.
        """
        wavelength, mags = self._read_csv(self.CSV_FILE, moon_phase)
        sed = self._calculate_photons_density(mags)
        sed = self._interpolate_spectral_distribution(
            wavelength, sed, object_wavelength)
        return sed
