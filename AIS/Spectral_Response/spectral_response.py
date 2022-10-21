"""
Spectral Respose Class
======================


The Spectral Response is an abstract class that represents the optical systems in the light path. 
Theses systems are the atmosphere, the telescope, and the SPARC4 instrument.
"""


import os
import pandas as pd
from scipy.interpolate import splev, splrep
import numpy as np
from numpy import ndarray


class Spectral_Response:
    """Spectral Response class

    The Spectral Response is an abstract class that represents the optical systems present in the light path 
    between the object and a detector on the ground.     

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """
    _BASE_PATH = os.path.join('AIS', 'Spectral_Response')
    _CSV_FILE_NAME = 'csv_file.csv'

    def __init__(self) -> None:
        """Initialize the class."""
        return

    @staticmethod
    def _read_csv_file(csv_file_name):
        ss = pd.read_csv(csv_file_name)
        return np.array(ss['Wavelength (nm)']), np.array(ss['Transmitance (%)'])

    @staticmethod
    def _interpolate_spectral_response(wavelength, spectral_response, obj_wavelength):
        spl = splrep(wavelength, spectral_response)
        spectral_response = splev(obj_wavelength, spl)
        return spectral_response

    def get_spectral_response(self, obj_wavelength: ndarray) -> ndarray:
        """Return the spectral response.

        Parameters
        ----------
        obj_wavelength: array like
            The wavelength interval, in nm, of the object.

        Yields
        ------
        spectral_response: array like
            The spectral response of the optical system.            
        """
        csv_file_name = os.path.join(self._BASE_PATH, self._CSV_FILE_NAME)
        sys_wavelength, spectral_response = self._read_csv_file(
            csv_file_name)
        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response, obj_wavelength)

        return spectral_response

    def apply_spectral_response(self, esd: ndarray, obj_wavelength: ndarray) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            esd: array like
                The Energy Spectral Distribution (ESD) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.

        Yields
        ------        
            reduced_esd: array like
                The Energy Spectral Distribution (ESD) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        spectral_response = self.get_spectral_response(obj_wavelength)
        reduced_esd = np.multiply(spectral_response, esd)

        return reduced_esd


class Telescope(Spectral_Response):
    """
    Telescope class

    This class extends the Spectral Response class, and it represents the spectral response of the
    1.6 m Perkin-Elmer telescope of the Picos dos Dias observatory.

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """

    _CSV_FILE_NAME = 'telescope.csv'

    def __init__(self) -> None:
        """Initialize the class."""
        super().__init__()
        return


class Atmosphere(Spectral_Response):
    """Atmosphere class

    This class extends the Spectral Reponse class, and it represents the atmosphere spectral response

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """

    _CSV_FILE_NAME = 'atmosphere.csv'

    def __init__(self) -> None:
        super().__init__()
        return

    @staticmethod
    def _read_csv_file(csv_file_name, sky_condition):
        ss = pd.read_csv(csv_file_name)
        return np.array(ss['Wavelength (nm)']), np.array(ss[sky_condition])/100

    def get_spectral_response(self, obj_wavelength: ndarray, air_mass: int | float, sky_condidition='photometric') -> ndarray:
        """Return the spectral response.

        Parameters
        ----------
        obj_wavelength: array like
            The wavelength interval, in nm, of the object.

        sky_condition: ['photometric', 'regular', 'good']
            The condition of the sky.
        air_mass: int, float
            The air mass.

        Yields
        ------
        spectral_response: array like
            The spectral response of the optical system.     

        """
        csv_file_name = os.path.join(self._BASE_PATH, self._CSV_FILE_NAME)
        sys_wavelength, extinction_coef = self._read_csv_file(
            csv_file_name, sky_condidition)
        spectral_response = 10 ** (-0.4 * extinction_coef * air_mass)
        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response, obj_wavelength)

        return spectral_response

    def apply_spectral_response(self, esd: ndarray, obj_wavelength: ndarray,
                                air_mass: int | float, sky_condition='photometric') -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            esd: array like
                The Energy Spectral Distribution (ESD) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.
            air_mass: int, float
                The air mass in the light path.
            sky_condition: ['photometric', 'regular', 'good']

        Yields
        ------        
            reduced_esd: array like
                The Energy Spectral Distribution (ESD) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        spectral_response = self.get_spectral_response(
            obj_wavelength, air_mass, sky_condition)
        reduced_esd = np.multiply(spectral_response, esd)

        return reduced_esd


class Channel(Spectral_Response):
    """
    Channel class

    This is an abstract class that extends the Spectral Response class. It represents one of the four
    SPARC4 channels.

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """

    def __init__(self) -> None:
        super().__init__()
        return
