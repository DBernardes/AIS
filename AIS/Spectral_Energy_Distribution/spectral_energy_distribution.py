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


class Spectral_Energy_Distribution:
    """Spectral Energy Distribution Class

    The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
    """

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


class Source(Spectral_Energy_Distribution):
    """Source Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the astronomical object.
    """
    EFFECT_WAVELENGTH = 550  # nm
    TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    S_0 = vega_fluxd.get()["Johnson V"].value*1e7  # W/(m.m2)
    NAME_SED_SPECTRAL_TYPE = {'O': 'uko5v.dat', 'B': 'ukb0i.dat', 'A': 'uka0i.dat',
                              'F': 'ukf0i.dat', 'G': 'ukg0i.dat', 'K': 'ukk0iii.dat', 'M': 'ukm0iii.dat'}

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
            spectral response and the wavelength of the object are obtained using a library of spectral types. The
            spectral types made available are: 'O', 'B', 'A', 'F', 'G', 'K', 'M'.  
            These spectrums are taken from the Library of Stellar Spectrum of the ESO, and they can be found at:
            https://www.eso.org/sci/observing/tools/standards/spectra/index.html.
            The level of the spectral response is adjusted using the effective flux. 
            The effective flux is calculated using the magnitude of the object in the V band.

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

        spectral_type : ['O', 'B', 'A', 'F', 'G', 'K', 'M'], optional
            The spectral type of the star that will be used to calculate the SED.
            This parameter is used only if the calculation_method is 'spectral_standard'.

        Returns
        -------
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
            wavelength, sed = self._read_spectral_library(spectral_type)
        else:
            raise ValueError(
                "The calculation_method must be 'blackbody' or 'spectral_library'.")

        effective_flux = self._calculate_effective_flux(
            magnitude)
        sed = sed * effective_flux
        return wavelength, sed

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature):
        wavelength *= 1e-9
        return 2 * pi * h * c ** 2 / wavelength ** 5 * 1 / (np.exp(h * c / (wavelength * k * temperature)) - 1)

    def _calculate_effective_flux(self, magnitude):
        return self.S_0*10**(-magnitude/2.5)*self.TELESCOPE_EFFECTIVE_AREA*self.EFFECT_WAVELENGTH*1e-9/h*c

    def _read_spectral_library(self, spectral_type):
        path = os.path.join(self.SPECTRAL_LIB_PATH,
                            self.NAME_SED_SPECTRAL_TYPE[spectral_type])
        sed = np.loadtxt(path)
        return sed[:, 0]/10, sed[:, 1]


class Sky(Spectral_Energy_Distribution):
    """Sky Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the emission of the sky. 
    """

    @staticmethod
    def _read_csv(file_path, value_name):
        ss = pd.read_csv(file_path)
        wavelenght = ss["wavelength"]
        value = ss[value_name]

        return wavelenght, value
