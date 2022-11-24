"""
Spectral Respose Class
======================


The Spectral Response is an abstract class that represents the optical systems in the light path. 
Theses systems are the atmosphere, the telescope, and the SPARC4 instrument.
"""


import os

import numpy as np
import pandas as pd
from numpy import ndarray
from scipy.interpolate import splev, splrep

__all__ = ['Atmosphere', 'Telescope', 'Channel']


class Spectral_Response:
    """Spectral Response class

    The Spectral Response is an abstract class that represents the optical systems present in the light path
    between the object and a detector on the ground.

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """

    _BASE_PATH = os.path.join("AIS", "Spectral_Response")
    _CSV_FILE_NAME = "csv_file.csv"

    def __init__(self) -> None:
        """Initialize the class."""
        return

    @staticmethod
    def _check_var_in_a_list(var, var_name, _list):
        if var not in _list:
            raise ValueError(
                f"The allowed values for the {var_name} are: {_list}")

    @staticmethod
    def _verify_var_in_interval(var, var_name, var_min=0, var_max=2 ** 32):
        if var <= var_min or var >= var_max:
            raise ValueError(
                f"The {var_name} must be in the interval [{var_min},{var_max}]: {var}."
            )

    @staticmethod
    def _read_csv_file(csv_file_name):
        ss = pd.read_csv(csv_file_name)
        return np.array(ss["Wavelength (nm)"]), np.array(ss["Transmitance (%)"]) / 100

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
        sys_wavelength, spectral_response = self._read_csv_file(csv_file_name)
        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response, obj_wavelength
        )

        return spectral_response

    def apply_spectral_response(self, sed: ndarray, obj_wavelength: ndarray) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            sed: array like
                The Spectral Energy Distribution (SED) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.

        Yields
        ------
            reduced_sed: array like
                The Spectral Energy Distribution (SED) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        spectral_response = self.get_spectral_response(obj_wavelength)
        reduced_sed = np.multiply(spectral_response, sed)

        return reduced_sed


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

    _CSV_FILE_NAME = "telescope.csv"

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

    _CSV_FILE_NAME = "atmosphere.csv"

    def __init__(self) -> None:
        super().__init__()
        return

    @staticmethod
    def _read_csv_file(csv_file_name, sky_condition):
        ss = pd.read_csv(csv_file_name)
        return np.array(ss["Wavelength (nm)"]), np.array(ss[sky_condition]) / 100

    def get_spectral_response(
        self,
        obj_wavelength: ndarray,
        air_mass: int | float,
        sky_condidition="photometric",
    ) -> ndarray:
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
            csv_file_name, sky_condidition
        )
        spectral_response = 10 ** (-0.4 * extinction_coef * air_mass)
        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response, obj_wavelength
        )

        return spectral_response

    def apply_spectral_response(
        self,
        sed: ndarray,
        obj_wavelength: ndarray,
        air_mass: int | float,
        sky_condition="photometric",
    ) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            sed: array like
                The Spectral Energy Distribution (SED) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.
            air_mass: int, float
                The air mass in the light path.
            sky_condition: ['photometric', 'regular', 'good']

        Yields
        ------
            reduced_sed: array like
                The Spectral Energy Distribution (SED) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        self._verify_var_in_interval(air_mass, "air_mass", var_max=3)
        self._check_var_in_a_list(sky_condition, "sky_condition", [
                                  "photometric", "regular", "good"])
        spectral_response = self.get_spectral_response(
            obj_wavelength, air_mass, sky_condition
        )
        reduced_sed = np.multiply(spectral_response, sed)

        return reduced_sed


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

    _POL_OPTICAL_COMPONENTS = {
        "polarizer": "polarizer.csv",
        "depolarizer": "depolarizer.csv",
        "retarder": "retarder.csv",
        "analyser": "analyser.csv",
    }

    _PHOT_OPTICAL_COMPONENTS = {
        "collimator": "collimator.csv",
        "dichroic": "dichroic.csv",
        "camera": "camera.csv",
        "ccd": "ccd.csv",
    }

    def __init__(self, channel_id: int | float) -> None:
        """Initialize the class.

        Parameters
        ----------

        channel_id: str
            The channel ID number.
        """

        super().__init__()
        self._channel_id = channel_id
        self._BASE_PATH = os.path.join(self._BASE_PATH, "channel")
        return

    def get_spectral_response(
        self, obj_wavelength: ndarray, csv_file_name: str
    ) -> ndarray:
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
        self._CSV_FILE_NAME = csv_file_name
        return super().get_spectral_response(obj_wavelength)

    def write_sparc4_operation_mode(self, acquisition_mode: str,
                                    calibration_wheel: str = 'empty',
                                    retarder_waveplate: str = 'half') -> None:
        """Write the operation mode of the SPARC4.

        Parameters
        ----------
        acquisition_mode: ['polarimetry', 'photometry']
            The acquisition mode of the channel. If the acquisition mode is polarimetry,
            the polarimetric configuration must be provided.

        calibration_wheel: ["polarizer", "depolarizer", "empty"], optional
            The position of the calibration wheel.            

        retarder: ["half", "quarter"], optional
            The waveplate for polarimetric measurements.
        """
        self.acquisition_mode = acquisition_mode
        self.calibration_wheel = calibration_wheel
        self.retarder_waveplate = retarder_waveplate
        self._verify_sparc4_operation_mode()
        return

    def apply_spectral_response(self, sed: ndarray, obj_wavelength: ndarray) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            sed: array like
                The Spectral Energy Distribution (SED) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.

        Yields
        ------
            reduced_sed: array like
                The Spectral Energy Distribution (SED) of the object subtracted from the loses
                related to the spectral response of the system.
        """

        if self.acquisition_mode == "polarimetry":
            sed = self._apply_polarimetric_spectral_response(
                sed, obj_wavelength)
        reduced_sed = self._apply_photometric_spectral_response(
            sed, obj_wavelength)

        return reduced_sed

    def _apply_photometric_spectral_response(self, sed: ndarray, obj_wavelength: ndarray) -> ndarray:
        for csv_file in self._PHOT_OPTICAL_COMPONENTS.values():
            if csv_file != "collimator.csv":
                csv_file = os.path.join(
                    f"Channel {self._channel_id}", csv_file)
            spectral_response = self.get_spectral_response(
                obj_wavelength, csv_file)
            sed = np.multiply(spectral_response, sed)

        return sed

    def _apply_polarimetric_spectral_response(self, sed: ndarray, obj_wavelength: ndarray) -> ndarray:
        _dict = {}
        cal_wheel = self.calibration_wheel
        if cal_wheel != "empty":
            _dict[cal_wheel] = self._POL_OPTICAL_COMPONENTS[cal_wheel]
        _dict['retarder'] = self._POL_OPTICAL_COMPONENTS['retarder']
        _dict['analyser'] = self._POL_OPTICAL_COMPONENTS['analyser']
        for csv_file in _dict.values():
            spectral_response = self.get_spectral_response(
                obj_wavelength, csv_file)
            sed = np.multiply(spectral_response, sed)

        return sed

    def _verify_sparc4_operation_mode(self) -> None:

        if self.acquisition_mode == "photometry":
            pass
        elif self.acquisition_mode == "polarimetry":
            self._check_var_in_a_list(
                self.calibration_wheel,
                "calibration wheel",
                [
                    "polarizer",
                    "depolarizer",
                    "empty",
                ],
            )

            self._check_var_in_a_list(
                self.retarder_waveplate,
                "retarder waveplate",
                ["half", "quarter"],
            )
        else:
            raise ValueError(
                f"The SPARC4 acquisition mode should be 'photometry' or 'polarimetry': {self.acquisition_mode}."
            )
        return


if __name__ == "__main__":
    pass
