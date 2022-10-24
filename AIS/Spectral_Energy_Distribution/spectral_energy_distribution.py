"""
Spectral Energy Distribution Class
==================================


The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
"""
from os import setgid
from numpy import ndarray
import numpy as np
from scipy.interpolate import splev, splrep
import pandas as pd
from scipy.constants import c, h, k


class Spectral_Energy_Distribution:
    """Spectral Energy Distribution Class

    The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
    """

    def __init__(self):
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


class Source(Spectral_Energy_Distribution):
    """Source Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the astronomical object.
    """

    TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    S_0 = 4e-2  # W/m2/m #duvida!

    def calculate_sed(self,
                      calculation_method: str,
                      wavelength_interval: tuple,
                      magnitude: int | float,
                      effect_wavelength: int | float,
                      temperature: int | float = 0,
                      spectral_class: str = '') -> ndarray:
        """Get the Spectral Energy Distribution of the astronomical object.

        Parameters
        ----------

        calculation_method : ['blackbody', 'spectral_standard']
            The method used to calculate the SED.
            If the user wants to use a blackbody SED, the calculation_method must be 'blackbody'.
            If the user wants to use a spectral standard SED, the calculation_method must be 'spectral_standard'.

            Explain the blackbody and spectral standard methods...

        wavelength_interval : tuple
            The wavelength interval, in nm, of the astronomical object.
            This parameter must be a tuple with three elements, where the first element is the initial wavelength,
            the second element is the final wavelength and the third element is the number of elements in the array

        magnitude : int | float
            The magnitude of the astronomical object in one of the BVRI bands.

        effect_wavelength : int | float
            The effective wavelength, in nm, of the band of the provided magnitude.


        temperature : int | float, optional
            The blackbody temperature of the astronomical object in Kelvin.

        spectral_class : str, optional
            The spectral class of the star that will be used to calculate the SED.

        Returns
        -------
            esd : ndarray
                The SED of the astronomical object in photons/s/m.
        """
        wavelength = np.linspace(
            wavelength_interval[0], wavelength_interval[1], wavelength_interval[2]) * 1e-9

        if calculation_method == 'blackbody':
            sed = self._calculate_sed_blackbody(wavelength, temperature)
        elif calculation_method == 'spectral_standard':  # duvida!
            sed = []
        else:
            raise ValueError(
                "The calculation_method must be 'user', 'blackbody' or 'spectral_standard'.")

        effective_flux = self._calculate_effective_flux(
            magnitude, effect_wavelength)
        sed = sed * effective_flux / max(sed)  # duvida!
        return sed

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature):
        # corrigir
        return 2 * h * c ** 2 / wavelength ** 5 * 1 / (np.exp(h * c / (wavelength * k * temperature)) - 1)

    def _calculate_effective_flux(self, magnitude, effective_wavelength):
        return self.S_0*10**(-magnitude/2.5)*self.TELESCOPE_EFFECTIVE_AREA*effective_wavelength*1e-9/h*c


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

    @staticmethod
    def _interpolate_spectral_distribution(wavelength, spectral_distribution, obj_wavelength):
        spl = splrep(wavelength, spectral_distribution)
        interpolated_spectral_distribution = splev(obj_wavelength, spl)
        return interpolated_spectral_distribution
