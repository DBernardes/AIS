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
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep, interp1d
from ._utils import POLARIZER_90_MATRIX, POLARIZER_MATRIX, calculate_retarder_matrix
from scipy.optimize import curve_fit

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
        spl = interp1d(wavelength, spectral_response,
                       bounds_error=False, fill_value=0, kind='cubic')

        return spl(obj_wavelength)

    def get_spectral_response(self, obj_wavelength: ndarray) -> ndarray:
        """Return the spectral response.

        Parameters
        ----------
        obj_wavelength: array like
            The wavelength interval, in nm, of the object.

        Returns
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
    _PRIMARY_MIRROR_ADJUSTMENT = 0.933
    _SECONDARY_MIRROR_ADJUSTMENT = 0.996

    def __init__(self) -> None:
        """Initialize the class."""
        super().__init__()
        return

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
        spectral_response = super().get_spectral_response(obj_wavelength)

        return spectral_response**2 * self._PRIMARY_MIRROR_ADJUSTMENT * self._SECONDARY_MIRROR_ADJUSTMENT


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
        return np.array(ss["Wavelength (nm)"]), np.array(ss[sky_condition])

    def get_spectral_response(
        self,
        obj_wavelength: ndarray,
        air_mass: int | float = 1,
        sky_condidition="photometric",
    ) -> ndarray:
        """Return the spectral response.

        Parameters
        ----------
        obj_wavelength: array like
            The wavelength interval, in nm, of the object.

            air_mass: int, float
            The air mass.

        sky_condition: ['photometric', 'regular', 'good']
            The condition of the sky.

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
            sys_wavelength, spectral_response, obj_wavelength)
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

    @staticmethod
    def _interpolate_spectral_response(wavelength, spectral_response, obj_wavelength):
        popt, _ = curve_fit(func, wavelength,
                            spectral_response)
        spectral_response = func(obj_wavelength, *popt)
        return spectral_response


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

    _POLARIZER_ANGLE = 0

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
            csv_file_name: str
                The name of the csv file that contains the spectral response of the optical component.

        Returns
        ------
            spectral_response: array like
                The spectral response of the optical system.
        """
        self._CSV_FILE_NAME = csv_file_name
        return super().get_spectral_response(obj_wavelength)

    def write_sparc4_operation_mode(self, acquisition_mode: str,
                                    calibration_wheel: str = '',
                                    retarder_waveplate: str = 'half',
                                    retarder_waveplate_angle: float = 0) -> None:
        """Write the operation mode of the SPARC4.

        Parameters
        ----------
        acquisition_mode: ['polarimetry', 'photometry']
            The acquisition mode of the channel. If the acquisition mode is polarimetry,
            the polarimetric configuration must be provided.

        calibration_wheel: ["polarizer", "depolarizer"], optional
            The position of the calibration wheel.

        retarder_waveplate: ["half", "quarter"], optional
            The waveplate for polarimetric measurements.

        retarder_waveplate_angle: float, optional
            The angle of the retarder waveplate in degrees.
            If the acquisition mode of SPARC4 is polarimetry, the retarder waveplate angle should be provided.
        """
        self.acquisition_mode = acquisition_mode
        self.calibration_wheel = calibration_wheel
        self.retarder_waveplate = retarder_waveplate
        self.retarder_waveplate_angle = retarder_waveplate_angle
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

        Returns
        ------
            reduced_sed: array like
                The Spectral Energy Distribution (SED) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        self.sed = sed
        self.obj_wavelength = obj_wavelength
        if self.acquisition_mode == "polarimetry":
            self._resize_sed()
            self._apply_polarimetric_spectral_response()
        self._apply_photometric_spectral_response()

        return self.sed

    def _apply_photometric_spectral_response(self) -> None:
        for csv_file in self._PHOT_OPTICAL_COMPONENTS.values():
            if csv_file != "collimator.csv":
                csv_file = os.path.join(
                    f"Channel {self._channel_id}", csv_file)
            spectral_response = self.get_spectral_response(
                self.obj_wavelength, csv_file)
            self.sed = np.multiply(spectral_response, self.sed)

        return

    def _apply_polarimetric_spectral_response(self) -> None:
        if self.calibration_wheel != '':
            self._apply_calibration_wheel()
        self._apply_retarder_waveplate()
        self._apply_analyser()
        return

    def _apply_calibration_wheel(self) -> None:
        spectral_response = self.get_spectral_response(
            self.obj_wavelength, self._POL_OPTICAL_COMPONENTS[self.calibration_wheel])
        self.sed[0, :] = np.multiply(spectral_response, self.sed[0])

        if self.calibration_wheel == "polarizer":
            POLARIZER_MATRIX = self._calc_polarizer_matrix(
                self._POLARIZER_ANGLE)
            self.sed = np.transpose([POLARIZER_MATRIX.dot(self.sed[:, i])
                                     for i in range(self.sed.shape[1])])
        elif self.calibration_wheel == "depolarizer":
            sed = self.sed
            self.sed = np.zeros((4, sed.shape[1]))
            self.sed[0, :] = sed[0, :]
        else:
            raise ValueError(
                f"The calibration wheel {self.calibration_wheel} is not valid.")

        return

    def _apply_retarder_waveplate(self) -> None:
        spectral_response = self.get_spectral_response(
            self.obj_wavelength, self._POL_OPTICAL_COMPONENTS['retarder'])
        self.sed[0, :] = np.multiply(spectral_response, self.sed[0])

        RETARDER_MATRIX = self._calc_retarder_matrix()
        self.sed = np.transpose([RETARDER_MATRIX.dot(self.sed[:, i])
                                 for i in range(self.sed.shape[1])])
        return

    def _apply_analyser(self) -> None:
        spectral_response = self.get_spectral_response(
            self.obj_wavelength, self._POL_OPTICAL_COMPONENTS['analyser'])
        self.sed[0, :] = np.multiply(spectral_response, self.sed[0])

        ORD_RAY_MATRIX = self._calc_polarizer_matrix(self._POLARIZER_ANGLE)
        temp_1 = np.transpose([ORD_RAY_MATRIX.dot(self.sed[:, i])
                               for i in range(self.sed.shape[1])])

        EXTRA_ORD_RAY_MATRIX = self._calc_polarizer_matrix(
            self._POLARIZER_ANGLE + 90)
        # temp_2 esta negativo !
        temp_2 = np.transpose([EXTRA_ORD_RAY_MATRIX.dot(self.sed[:, i])
                               for i in range(self.sed.shape[1])])

        self.sed = np.stack((temp_1[0], temp_2[0]))

    def _apply_polarimetric_spectral_response_1(self, sed: ndarray, obj_wavelength: ndarray) -> ndarray:
        _dict = {}
        cal_wheel = self.calibration_wheel
        if cal_wheel != '':
            _dict[cal_wheel] = self._POL_OPTICAL_COMPONENTS[cal_wheel]
        _dict['retarder'] = self._POL_OPTICAL_COMPONENTS['retarder']
        _dict['analyser'] = self._POL_OPTICAL_COMPONENTS['analyser']

        for csv_file in _dict.values():
            spectral_response = self.get_spectral_response(
                obj_wavelength, csv_file)
            sed[0, :] = np.multiply(spectral_response, sed[0])

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
                    "",
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

    def _resize_sed(self):
        sed = np.array(self.sed)
        if sed.ndim == 1:
            n = len(sed)
            temp = sed
            sed = np.zeros((4, n))
            sed[0, :] = temp
        self.sed = sed
        return

    @staticmethod
    def _calc_polarizer_matrix(polarizer_angle: float) -> ndarray:
        pol_angle = np.radians(polarizer_angle)
        POLARIZER_MATRIX = 0.5 * np.asarray(
            [
                [1, cos(2 * pol_angle), sin(2 * pol_angle), 0],
                [
                    cos(2 * pol_angle),
                    cos(2 * pol_angle) ** 2,
                    cos(2 * pol_angle) * sin(2 * pol_angle),
                    0,
                ],
                [
                    sin(2 * pol_angle),
                    cos(2 * pol_angle) * sin(2 * pol_angle),
                    sin(2 * pol_angle) ** 2,
                    0,
                ],
                [0, 0, 0, 0],
            ]
        )
        return POLARIZER_MATRIX

    def _calc_retarder_matrix(self):
        if self.retarder_waveplate == "half":
            phase_difference = np.radians(90)
        elif self.retarder_waveplate == "quarter":
            phase_difference = np.radians(45)
        else:
            raise ValueError(
                f"The retarder waveplate {self.retarder_waveplate} is not valid."
            )
        ret_angle = np.radians(self.retarder_waveplate_angle)
        RETARDER_MATRIX = np.asarray(
            [
                [1, 0, 0, 0],
                [
                    0,
                    cos(2 * ret_angle) ** 2
                    + sin(2 * ret_angle) ** 2 * cos(phase_difference),
                    cos(2 * ret_angle)
                    * sin(2 * ret_angle)
                    * (1 - cos(phase_difference)),
                    -sin(2 * ret_angle) * sin(phase_difference),
                ],
                [
                    0,
                    cos(2 * ret_angle)
                    + sin(2 * ret_angle) * (1 - cos(phase_difference)),
                    sin(2 * ret_angle) ** 2
                    + cos(2 * ret_angle) ** 2 * cos(phase_difference),
                    cos(2 * ret_angle) * sin(phase_difference),
                ],
                [
                    0,
                    sin(2 * ret_angle) * sin(phase_difference),
                    -cos(2 * ret_angle) * sin(phase_difference),
                    cos(phase_difference),
                ],
            ]
        )

        return RETARDER_MATRIX


def func(x, a, b, c, d):
    return a*x**3 + b*x**2 + c*x + d
