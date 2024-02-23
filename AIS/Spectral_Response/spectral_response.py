"""
Spectral Respose Class
======================


The Spectral Response is an abstract class that represents the optical systems in the light path.
Theses systems are the atmosphere, the telescope, and the SPARC4 instrument.
"""

import os
from copy import copy
from math import atan2, sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import ndarray
from scipy.interpolate import (
    PchipInterpolator,
    UnivariateSpline,
    interp1d,
    splev,
    splrep,
)
from scipy.optimize import curve_fit

from AIS.Spectral_Response._utils import (
    apply_matrix,
    calculate_depolarizer_matrix,
    calculate_polarizer_matrix,
    calculate_retarder_matrix,
)

__all__ = ["Atmosphere", "Telescope", "Channel"]


class Spectral_Response:
    """Spectral Response class

    The Spectral Response is an abstract class that represents the optical systems present in the light path
    between the object and a detector on the ground.

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """

    BASE_PATH = os.path.join("AIS", "Spectral_Response")
    CSV_FILE_NAME = "csv_file.csv"

    def __init__(self) -> None:
        """Initialize the class."""
        return

    @staticmethod
    def _check_var_in_a_list(var, var_name, _list) -> None:
        if var not in _list:
            raise ValueError(f"The allowed values for the {var_name} are: {_list}")

    @staticmethod
    def _verify_var_in_interval(var, var_name, var_min=0, var_max=2**32) -> None:
        if var <= var_min or var >= var_max:
            raise ValueError(
                f"The {var_name} must be in the interval [{var_min},{var_max}]: {var}."
            )

    @staticmethod
    def _read_csv_file(csv_file_name) -> tuple[ndarray]:
        ss = pd.read_csv(csv_file_name, dtype=np.float64)
        return np.array(ss["Wavelength (nm)"]), np.array(ss["Transmitance (%)"]) / 100

    def _interpolate_spectral_response(self, wavelength, spectral_response) -> ndarray:
        b = PchipInterpolator(wavelength, spectral_response)
        spectral_response = b(self.obj_wavelength)

        return spectral_response

    def get_spectral_response(self, obj_wavelength: ndarray) -> ndarray:
        """Return the spectral response.

        Parameters
        ----------
        obj_wavelength: array like.
            Wavelength interval of the object in nm.

        Returns
        ------
        spectral_response: array like
            The spectral response of the optical system.
        """
        self.obj_wavelength = obj_wavelength
        csv_file_name = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)
        sys_wavelength, spectral_response = self._read_csv_file(csv_file_name)
        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response
        )

        return spectral_response

    def apply_spectral_response(self, obj_wavelength: ndarray, sed: ndarray) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
        obj_wavelength: array like
            The wavelength interval, in nm, of the object.
        sed: array like in the (4, n) format
            The Spectral Energy Distribution (SED) of the object.

        Yields
        ------
        reduced_sed: array like
            The Spectral Energy Distribution (SED) of the object subtracted from the loses
            related to the spectral response of the system.
        """
        spectral_response = self.get_spectral_response(obj_wavelength)
        sed = sed * spectral_response

        return sed


class Telescope(Spectral_Response):
    CSV_FILE_NAME = "telescope.csv"
    PRIMARY_MIRROR_ADJUSTMENT = 0.933
    SECONDARY_MIRROR_ADJUSTMENT = 0.996

    def __init__(self) -> None:
        """
        Telescope class

        This class extends the Spectral Response class, and it represents the spectral response of the
        1.6 m Perkin-Elmer telescope of the Picos dos Dias observatory.

        Yields
        ------
            spectral_response: array_like
                The spectral response of the optical system.
        """
        super().__init__()
        return

    def get_spectral_response(self, obj_wavelength) -> ndarray:
        """Return the spectral response.

        Parameters
        ---------
        obj_wavelength: array like
            Wavelength interval of the object in nm.

        Returns
        ------

        spectral_response: array like
            The spectral response of the optical system.
        """
        spectral_response = super().get_spectral_response(obj_wavelength)

        return (
            spectral_response**2
            * self.PRIMARY_MIRROR_ADJUSTMENT
            * self.SECONDARY_MIRROR_ADJUSTMENT
        )


class Atmosphere(Spectral_Response):
    CSV_FILE_NAME = "atmosphere_profile.csv"
    ATM_EXTINCTION_FILE = "willton.csv"

    def __init__(self) -> None:
        """Atmosphere class

        This class extends the Spectral Reponse class, and it represents the atmosphere spectral response

        Yields
        ------
            spectral_response: array_like
                The spectral response of the optical system.
        """
        super().__init__()
        return

    def get_spectral_response(
        self,
        obj_wavelength: ndarray,
        air_mass: int | float = 1,
        sky_condition="photometric",
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
        ss = pd.read_csv(
            os.path.join(self.BASE_PATH, self.ATM_EXTINCTION_FILE), dtype=np.float64
        )
        spectral_response = 10 ** (-0.4 * air_mass * ss[sky_condition])
        spectral_response = self._interpolate_spectral_response(
            ss["Wavelength (nm)"], spectral_response, obj_wavelength
        )
        spectral_response[spectral_response > 1] = 1

        return spectral_response

    def get_spectral_response_2(
        self,
        obj_wavelength: ndarray,
        air_mass: int | float = 1,
        sky_condition="photometric",
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
        ss = pd.read_csv(
            os.path.join(self.BASE_PATH, self.ATM_EXTINCTION_FILE), dtype=np.float64
        )
        extinction_coef = 10 ** (-0.4 * air_mass * ss[sky_condition])
        popt_M1, _ = curve_fit(self._func, ss["Wavelength (nm)"], extinction_coef)

        return super().get_spectral_response(obj_wavelength) * popt_M1[0]

    def apply_spectral_response(
        self,
        obj_wavelength: ndarray,
        sed: ndarray,
        air_mass: int | float,
        sky_condition="photometric",
    ) -> ndarray:
        """Apply the spectral response.

        Parameters
        ----------
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.
            sed: array like in the (4, n) format
                The Spectral Energy Distribution (SED) of the object.
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
        self._check_var_in_a_list(
            sky_condition, "sky_condition", ["photometric", "regular", "good"]
        )
        spectral_response = self.get_spectral_response(
            obj_wavelength, air_mass, sky_condition
        )

        sed = spectral_response * sed

        return sed

    def _func(self, wavelength_interv: ndarray, C: float) -> ndarray:
        """Profile for the atmosphere spectral response

        Parameters
        ----------
        wavelength_interv : ndarray
            wavelength interval
        C : float
            Adjustment constant

        Returns
        -------
        ndarray
            Adjusted profile
        """
        csv_file_name = os.path.join(self.BASE_PATH, self.CSV_FILE_NAME)
        sys_wavelength, spectral_response = self._read_csv_file(csv_file_name)
        spl = interp1d(
            sys_wavelength, spectral_response, bounds_error=False, kind="cubic"
        )
        return spl(wavelength_interv) * C

    @staticmethod
    def _interpolate_spectral_response(
        wavelength, spectral_response, obj_wavelength
    ) -> ndarray:
        spl = interp1d(
            wavelength,
            spectral_response,
            bounds_error=False,
            fill_value="extrapolate",
            kind="linear",
        )

        return spl(obj_wavelength)


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

    POL_OPTICAL_COMPONENTS = {
        "analyzer": "analyzer.csv",
        "retarder": "retarder.csv",
    }
    PHOT_OPTICAL_COMPONENTS = {
        "collimator": "collimator.csv",
        "dichroic": "dichroic.csv",
        "camera": "camera.csv",
        "ccd": "ccd.csv",
    }

    POLARIZER_ANGLE = 0

    def __init__(self, channel_id: int) -> None:
        """Initialize the class.

        Parameters
        ----------

        channel_id: int
            The channel ID number.
        """

        super().__init__()
        self.channel_id = channel_id
        self.BASE_PATH = os.path.join(self.BASE_PATH, "channel")
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
        self.CSV_FILE_NAME = csv_file_name
        return super().get_spectral_response(obj_wavelength)

    def write_sparc4_operation_mode(
        self,
        acquisition_mode: str,
        calibration_wheel: str = "",
        retarder_waveplate: str = "half",
        retarder_waveplate_angle: float = 0,
    ) -> None:
        """Write the operation mode of the SPARC4.

        Parameters
        ----------
        acquisition_mode: ['polarimetry', 'photometry']
            The acquisition mode of the channel. If the acquisition mode is polarimetry,
            the polarimetric configuration must be provided.

        calibration_wheel: ["polarizer", "ideal-polarizer", "ideal-depolarizer", "depolarizer"], optional
            The optical component of the calibration wheel.
            This parameter provides to the user the options of using the real or the ideal versions
            of the optical components of the calibration wheel. Based on the provided value, AIS
            will apply the correspondent Stoke matrix. It should be highlighted that the transmission
            of the optical component still be applied.

        retarder_waveplate: ["ideal-half", "half", "ideal-quarter", "quarter"], optional
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
            sed: array like in the (4, n) format
                The Spectral Energy Distribution (SED) of the object.
            obj_wavelength: array like
                The wavelength interval, in nm, of the object.

        Returns
        ------
            reduced_sed: array like
                The Spectral Energy Distribution (SED) of the object subtracted from the loses
                related to the spectral response of the system.
        """
        self.sed = sed[0]
        self.obj_wavelength = obj_wavelength
        if self.acquisition_mode == "polarimetry":
            self.sed = sed
            self._apply_polarimetric_spectral_response()

        self._apply_photometric_spectral_response()

        return self.sed

    def _apply_photometric_spectral_response(self) -> None:
        for csv_file in self.PHOT_OPTICAL_COMPONENTS.values():
            if csv_file != "collimator.csv":
                csv_file = os.path.join(f"Channel {self.channel_id}", csv_file)
            spectral_response = self.get_spectral_response(
                self.obj_wavelength, csv_file
            )
            self.sed = spectral_response * self.sed
        #     plt.plot(self.obj_wavelength, self.sed[0], label=csv_file)
        #     plt.plot(self.obj_wavelength, self.sed[1], label=csv_file)
        #     ord = np.trapz(self.sed[0], self.obj_wavelength * 1e-9)
        #     extra = np.trapz(self.sed[1], self.obj_wavelength * 1e-9)
        #     print((ord - extra) / (ord + extra))
        # plt.show()
        return

    def _apply_polarimetric_spectral_response(self) -> None:
        for csv_file in self.POL_OPTICAL_COMPONENTS.values():
            spectral_response = self.get_spectral_response(
                self.obj_wavelength, csv_file
            )
            self.sed = spectral_response * self.sed
        #  ? Pode isso: fiz o teste uma vez e deu a mesma coisa

        if self.calibration_wheel != "":
            self._apply_calibration_wheel()
        self._apply_retarder_waveplate()
        self._apply_analyzer()

        return

    def _apply_calibration_wheel(self) -> None:
        if "depolarizer" in self.calibration_wheel:
            self._apply_depolarizer()
        elif "polarizer" in self.calibration_wheel:
            self._apply_polarizer()
        else:
            raise ValueError(
                f"The calibration wheel {self.calibration_wheel} is not valid."
            )

        return

    def _apply_polarizer(self) -> None:
        spectral_response = self.get_spectral_response(
            self.obj_wavelength, "polarizer.csv"
        )
        if self.calibration_wheel == "ideal-polarizer":
            self._apply_ideal_polarizer(spectral_response)
        elif self.calibration_wheel == "polarizer":
            self._apply_real_polarizer(spectral_response)
        else:
            raise ValueError(
                f"A wrong value was provided for the polarizer: {self.calibration_wheel}"
            )

    def _apply_ideal_polarizer(self, spectral_response) -> None:
        self.sed = spectral_response * self.sed
        polarizer_matrix = calculate_polarizer_matrix(self.POLARIZER_ANGLE)
        self.sed = apply_matrix(polarizer_matrix, self.sed)

    def _apply_real_polarizer(self, spectral_response) -> None:
        contrast_ratio = self._get_spectral_response_custom(
            "polarizer_contrast_ratio.csv", "Contrast ratio"
        )
        for idx, transx in enumerate(spectral_response):
            contrast = contrast_ratio[idx]
            transy = transx / contrast
            total_transmission = sqrt(transx + transy)
            theta = np.rad2deg(atan2(1, sqrt(contrast)))
            polarizer_matrix = total_transmission * calculate_polarizer_matrix(
                self.POLARIZER_ANGLE + theta
            )
            self.sed[:, idx] = np.transpose(polarizer_matrix.dot(self.sed[:, idx]))

    def _apply_real_polarizer_1(self, spectral_response) -> None:
        contrast_ratio = self._get_spectral_response_custom(
            "polarizer_contrast_ratio.csv", "Contrast ratio"
        )
        for idx, transmission in enumerate(spectral_response):
            contrast = contrast_ratio[idx]
            theta = np.rad2deg(atan(1 / sqrt(contrast)))
            total_transmission = transmission * (1 + 1 / contrast)
            polarizer_matrix = calculate_polarizer_matrix(self.POLARIZER_ANGLE + theta)
            self.sed[:, idx] = total_transmission * np.transpose(
                polarizer_matrix.dot(self.sed[:, idx])
            )
            self.sed[1, idx] *= 1 - 1 / contrast  # ? ta certo isso

    def _apply_depolarizer(self) -> None:
        spectral_response = self.get_spectral_response(
            self.obj_wavelength, "depolarizer.csv"
        )
        if self.calibration_wheel == "ideal-depolarizer":
            self._apply_ideal_depolarizer(spectral_response)
        elif self.calibration_wheel == "depolarizer":
            self._apply_real_depolarizer(spectral_response)
        else:
            raise ValueError(
                f"A wrong value was provided for the depolarizer: {self.calibration_wheel}"
            )

        return

    def _apply_ideal_depolarizer(self, spectral_response) -> None:
        tmp = self.sed
        self.sed = np.zeros((4, tmp.shape[1]))
        self.sed[0] = tmp[0] * spectral_response

    def _apply_real_depolarizer(self, spectral_response) -> None:
        self.sed = spectral_response * self.sed
        for idx, wavelength in enumerate(self.obj_wavelength):
            depolarizer_matrix = calculate_depolarizer_matrix(wavelength)
            self.sed[:, idx] = np.transpose(depolarizer_matrix.dot(self.sed[:, idx]))

    def _apply_retarder_waveplate(self) -> None:
        if self.retarder_waveplate[:5] == "ideal":
            self._apply_ideal_waveplate()
        elif self.retarder_waveplate in ["half", "quarter"]:
            self._apply_real_waveplate()
        else:
            raise ValueError(
                f"A wrong value was provided for the depolarizer: {self.retarder_waveplate}"
            )
        return

    def _apply_real_waveplate(self) -> None:
        phase_difference = (
            self._get_spectral_response_custom(
                f"retarder_phase_diff_{self.retarder_waveplate}.csv", "Retardance"
            )
            * 360
        )
        # phase_difference[np.where(phase_difference > 1)] = 89
        # plt.plot(self.obj_wavelength, phase_difference)
        # plt.show()
        for idx, phase in enumerate(phase_difference):
            RETARDER_MATRIX = calculate_retarder_matrix(
                phase, self.retarder_waveplate_angle
            )
            self.sed[:, idx] = np.transpose(RETARDER_MATRIX.dot(self.sed[:, idx]))

    def _apply_ideal_waveplate(self) -> None:
        if "half" in self.retarder_waveplate:
            phase = 180
        elif "quarter" in self.retarder_waveplate:
            phase = 90
        else:
            raise ValueError(
                f"A wrong value has been provided for the retarder waveplate: {self.retarder_waveplate}"
            )
        RETARDER_MATRIX = calculate_retarder_matrix(
            phase, self.retarder_waveplate_angle
        )
        self.sed = apply_matrix(RETARDER_MATRIX, self.sed)

    def _apply_analyzer(self) -> None:
        ORD_RAY_MATRIX = calculate_polarizer_matrix(self.POLARIZER_ANGLE)
        temp_1 = apply_matrix(ORD_RAY_MATRIX, copy(self.sed))
        EXTRA_ORD_RAY_MATRIX = calculate_polarizer_matrix(self.POLARIZER_ANGLE + 90)
        temp_2 = apply_matrix(EXTRA_ORD_RAY_MATRIX, copy(self.sed))
        self.sed = np.stack((temp_1[0], temp_2[0]))

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
                    "ideal-polarizer",
                    "ideal-depolarizer",
                    "",
                ],
            )

            self._check_var_in_a_list(
                self.retarder_waveplate,
                "retarder waveplate",
                ["half", "quarter", "ideal-half", "ideal-quarter"],
            )
        else:
            raise ValueError(
                f"The SPARC4 acquisition mode should be 'photometry' or 'polarimetry': {self.acquisition_mode}."
            )
        return

    def _get_spectral_response_custom(self, csv_file: str, column_name: str) -> ndarray:
        csv_file_name = os.path.join(self.BASE_PATH, csv_file)
        ss = pd.read_csv(csv_file_name)
        sys_wavelength, spectral_response = np.asarray(
            ss["Wavelength (nm)"]
        ), np.asarray(ss[column_name])

        spectral_response = self._interpolate_spectral_response(
            sys_wavelength, spectral_response
        )

        return spectral_response

    def _interpolate_spectral_response(self, wavelength, spectral_response) -> ndarray:
        spl = splrep(
            wavelength,
            spectral_response,
        )
        spectral_response = splev(self.obj_wavelength, spl)
        spectral_response[spectral_response < 0] = 0
        return spectral_response
