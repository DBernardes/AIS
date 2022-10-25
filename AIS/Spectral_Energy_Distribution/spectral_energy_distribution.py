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


class Spectral_Energy_Distribution:
    """Spectral Energy Distribution Class

    The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
    """

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
    EFFECT_WAVELENGTH = 550e-9
    TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    S_0 = vega_fluxd.get()["Johnson V"].value*1e7  # W/(m.m2)

    def calculate_sed(self,
                      calculation_method: str,
                      magnitude: int | float,
                      wavelength_interval: tuple = (),
                      temperature: int | float = 0,
                      spectral_type: str = '') -> ndarray:
        """Get the Spectral Energy Distribution of the astronomical object.

        Parameters
        ----------

        calculation_method : ['blackbody', 'spectral_standard']
            The method used to calculate the SED.
            If the user wants to use a blackbody SED, the calculation_method must be 'blackbody'.
            If the user wants to use a spectral standard SED, the calculation_method must be 'spectral_standard'.

            Explain the blackbody and spectral standard methods...


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
                The SED of the astronomical object in photons/s/m.
        """

        wv = wavelength_interval
        if calculation_method == 'blackbody':
            wavelength = np.linspace(wv[0], wv[1], wv[2]) * 1e-9
            sed = self._calculate_sed_blackbody(wavelength, temperature)
        elif calculation_method == 'spectral_standard':  # duvida!
            wavelength = np.linspace(350, 1100, 100)
            sed = np.linspace(100, 1000, 100)
        else:
            raise ValueError(
                "The calculation_method must be 'blackbody' or 'spectral_standard'.")

        normalization_flux = self._interpolate_spectral_distribution(
            wavelength, sed, self.EFFECT_WAVELENGTH)
        effective_flux = self._calculate_effective_flux(
            magnitude)
        sed = sed * effective_flux / normalization_flux
        return sed

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature):
        # corrigir
        return 2 * h * c ** 2 / wavelength ** 5 * 1 / (np.exp(h * c / (wavelength * k * temperature)) - 1)

    def _calculate_effective_flux(self, magnitude):
        return self.S_0*10**(-magnitude/2.5)*self.TELESCOPE_EFFECTIVE_AREA*self.EFFECT_WAVELENGTH/h*c


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
