"""
Spectral Energy Distribution Class
==================================


The Spectral Energy Distribtution is an abstract class that represents the sky and the 
source classes.
"""

from numpy import ndarray
import numpy as np
from scipy.interpolate import splev, splrep, interp1d
import pandas as pd
from scipy.constants import c, h, k
from sbpy.calib import vega_fluxd
from math import pi
import os
from sys import exit
from ..Spectral_Response._utils import calculate_polarizer_matrix, apply_matrix
from copy import copy
__all__ = ['Source', 'Sky']


class Spectral_Energy_Distribution:
    """Spectral Energy Distribution Class

    The Spectral Energy Distribtution is an abstract class that represents the sky and the source classes.
    """
    _EFFECT_WAVELENGTH = 555.6  # nm
    _TELESCOPE_EFFECTIVE_AREA = 0.804  # m2

    # https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
    # https://cass.ucsd.edu/archive/physics/ph162/mags.html
    _S_0 = 3.658e-2  # W/(m.m2)
    _BASE_PATH = os.path.join(
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
    def _interpolate_spectral_distribution(wavelength, spectral_response, obj_wavelength):
        # spl = interp1d(wavelength, spectral_response,
        #               bounds_error=False, fill_value='extrapolate', kind='cubic')
        spl = splrep(wavelength, spectral_response)
        interpolated_spectral_distribution = splev(obj_wavelength, spl)
        return interpolated_spectral_distribution  # spl(obj_wavelength)

    def _calculate_photons_density(self, magnitude) -> float:
        return self._S_0*10**(-magnitude/2.5)*self._TELESCOPE_EFFECTIVE_AREA*self._EFFECT_WAVELENGTH*1e-9/(h*c)


class Source(Spectral_Energy_Distribution):
    """Source Class

    This class inherits from the Spectral_Energy_Distribution class, and it represents the astronomical object.
    
    Example
    -------
    
    src = Source()
    
    wv, sed = src.calculate_sed(calculation_method='blackbody', 
                            magnitude=12, 
                            wavelength_interval=(400, 1100, 100), 
                            temperature=5700)
    """

    def __init__(self):
        self.SPECTRAL_LIB_PATH = os.path.join(
            self._BASE_PATH, 'Spectral_Library')
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
                wavelength, sed, self._EFFECT_WAVELENGTH)
            sed /= normalization_flux
        elif calculation_method == 'spectral_library':
            wavelength, sed = self._read_spectral_library(
                spectral_type.lower())
        else:
            raise ValueError(
                "The calculation_method must be 'blackbody' or 'spectral_library'.")

        effective_flux = self._calculate_photons_density(
            magnitude)

        n = len(sed)
        self.sed = np.zeros((4, n))
        self.sed[0] = sed * effective_flux
        self.wavelength = wavelength

        return self.wavelength, self.sed


    def apply_polarization(self, polarization_mode:str ='linear', pol_angle:float=0, percent_pol:float=100) -> ndarray:
        """Apply polarization to SED.

        Parameters
        ----------
            polarization_mode: ['linear', 'circular'], optional
                Polarization mode.
            pol_angle: float, optional
                Polarization angle in degrees. If the selected polarization mode were linear, the polarization angle must be provided.
            percent_pol: float, optional
                Percentage of polarization.

        Returns
        -------
            sed: ndarray
                Polarized SED.
        """
        percent_pol /= 100       
        if polarization_mode == 'linear':
            polarizer_matrix = calculate_polarizer_matrix(pol_angle)
            polarized_sed = apply_matrix(polarizer_matrix, copy(self.sed)) * 2
            self.sed = percent_pol * polarized_sed + (1 - percent_pol) * self.sed
        else:
            raise ValueError(f'Unknow polarization mode: {polarization_mode}')
        
        return self.sed

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature):
        return 2 * pi * h * c**2 / (wavelength*1e-9)**5 * 1 / (np.exp(h*c/(wavelength * 1e-9 * k * temperature))-1)

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
        file_name = os.path.join(self._BASE_PATH, file_name)
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
        temp = self._interpolate_spectral_distribution(
            wavelength, sed, object_wavelength)
        sed = np.zeros((4, object_wavelength.shape[0]))
        sed[0] = temp
        return sed
