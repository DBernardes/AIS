"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""


import os
from cmath import phase

import numpy as np
import pandas as pd
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep

# from sys import exit


class Abstract_SPARC4_Spectral_Response:

    """Abstract class of the SPARC4 spectral response."""

    _CHANNEL_ID = 0
    _DIR_PATH = os.path.join("AIS", "SPARC4_Spectral_Response")
    _THETA_POL = np.deg2rad(0)
    _POLARIZER_MATRIX = 0.5 * np.asarray(
        [
            [1, cos(2 * _THETA_POL), sin(2 * _THETA_POL), 0],
            [
                cos(2 * _THETA_POL),
                cos(2 * _THETA_POL) ** 2,
                cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
                0,
            ],
            [
                sin(2 * _THETA_POL),
                cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
                sin(2 * _THETA_POL) ** 2,
                0,
            ],
            [0, 0, 0, 0],
        ]
    )

    _THETA_POL = np.deg2rad(_THETA_POL + pi / 2)
    _POLARIZER_90_MATRIX = 0.5 * np.asarray(
        [
            [1, cos(2 * _THETA_POL), sin(2 * _THETA_POL), 0],
            [
                cos(2 * _THETA_POL),
                cos(2 * _THETA_POL) ** 2,
                cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
                0,
            ],
            [
                sin(2 * _THETA_POL),
                cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
                sin(2 * _THETA_POL) ** 2,
                0,
            ],
            [0, 0, 0, 0],
        ]
    )

    def __init__(self, wavelength_interval):
        """
        Initialize the class.

        Parameters
        ----------

        wavelength_interval: array-like.
            Wavelength interval of the specific flux.
        """

        self.wavelength_interval_len = len(wavelength_interval)
        self.wavelength_interval = wavelength_interval

    def get_channel_ID(self):
        """
        Return the chanel ID.

        Returns
        -------
        _CHANNEL_ID : [1, 2, 3, 4]
            Channel ID
        """
        return self._CHANNEL_ID

    def write_specific_photons_per_second(self, specific_photons_per_second):
        """
        Write the specific photons per second.

        This function writes the specific number of photons per second of the object in the class.

        Parameters
        ----------

        specific_photons_per_second: array-like
            Photons per second of the object, as a function of the wavelength.
        """
        self.specific_photons_per_second = specific_photons_per_second

    def read_specific_photons_per_second(self):
        """Read the specific photons per second array."""
        return self.specific_photons_per_second

    def apply_calibration_wheel(self, component):
        """Apply calibration wheel spectral response and polarimetric manipulation.

        Parameters
        ----------

        component: [polarizer, depolarizer, empty]
            The component of the calibration wheel.
        """
        if component == "polarizer":
            self._apply_polarizer()
        elif component == "depolarizer":
            self._apply_depolarizer()
        elif component == "empty":
            pass

    def _apply_polarizer(self):
        self._multiply_matrices(self._POLARIZER_MATRIX)
        file = os.path.join(self._DIR_PATH, "polarizer" + ".csv")
        self._apply_optical_component_transmission(file)

    def _apply_depolarizer(self):
        temp = self.specific_photons_per_second
        self.specific_photons_per_second = [np.zeros((4, self.wavelength_interval_len))]
        self.specific_photons_per_second[0][0] = temp[0][0]
        file = os.path.join(self._DIR_PATH, "depolarizer" + ".csv")
        self._apply_optical_component_transmission(file)

    def apply_retarder(self, type):
        """Apply retarder spectral response and polarimetric manipulation.

        Parameters
        ----------

        type: [half, quarter]
            The type of the retarder used to manipulate the incoming light
        """
        file = os.path.join(self._DIR_PATH, "retarder_phase_diff_" + type + ".csv")
        ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
        retarder_wavelength = ss["(nm)"]
        retardance = ss["(waves)"]
        retardance = self._calculate_spline(retardance, retarder_wavelength)
        retardance = [val * 360 for val in retardance]

        for idx, value in enumerate(retardance):
            retarder_matrix = self._calculate_retarder_matrix(value)
            self.specific_photons_per_second[0][:, idx] = np.dot(
                retarder_matrix, self.specific_photons_per_second[0][:, idx]
            )
        file = os.path.join(self._DIR_PATH, "retarder_transmitance.csv")
        self._apply_optical_component_transmission(file)

    def apply_analyser(self):
        """Apply analyser spectral response and polarimetric manipulation"""
        temp = self.specific_photons_per_second
        file = os.path.join(self._DIR_PATH, "analyser.csv")
        self._multiply_matrices(self._POLARIZER_MATRIX)
        self._apply_optical_component_transmission(file)
        ordinary_ray = self.specific_photons_per_second[0]

        self.specific_photons_per_second = temp
        self._multiply_matrices(self._POLARIZER_90_MATRIX)
        self._apply_optical_component_transmission(file)
        extra_ordinary_ray = self.specific_photons_per_second[0]

        self.specific_photons_per_second = [ordinary_ray, extra_ordinary_ray]

    def apply_collimator(self):
        """Collimator spectral response."""
        file = os.path.join(self._DIR_PATH, "collimator.csv")
        self._apply_optical_component_transmission(file)

    def apply_dichroic(self):
        """
        Apply the dichroic spectral response.


        This functions applies the spectral response of the two
        dichroics that compose each channel.
        """
        for i in [1, 2]:
            file = os.path.join(
                self._DIR_PATH, f"Channel {self._CHANNEL_ID}", f"dichroic_{i}.csv"
            )
            self._apply_optical_component_transmission(file)

    def apply_camera(self):
        """Apply the camera spectral response."""
        file = os.path.join(self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "camera.csv")
        self._apply_optical_component_transmission(file)

    def apply_ccd(self):
        """Apply the ccd spectral response."""
        file = os.path.join(self._DIR_PATH, f"Channel {self._CHANNEL_ID}", "ccd.csv")
        self._apply_optical_component_transmission(file)

    def _read_spreadsheet(self, file):
        ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
        wavelength = ss["(nm)"]
        transmitance = ss["(%)"] / 100
        return wavelength, transmitance

    def _multiply_matrices(self, matrix):
        temp = []
        for array in self.specific_photons_per_second:
            for idx, value in enumerate(array[0]):
                array[:, idx] = np.dot(matrix, array[:, idx])
            temp.append(array)
        self.specific_photons_per_second = temp

    def _calculate_spline(self, value, component_wavelength_interv):
        spl = splrep(component_wavelength_interv, value)
        value = splev(self.wavelength_interval, spl)
        return value

    def _apply_optical_component_transmission(self, file):
        component_wavelength_interv, component_transmitance = self._read_spreadsheet(
            file
        )
        transmitance = self._calculate_spline(
            component_transmitance, component_wavelength_interv
        )
        temp = []
        for array in self.specific_photons_per_second:
            array[0] = np.multiply(array[0], transmitance)
            temp.append(array)
        self.specific_photons_per_second = temp

    def _calculate_retarder_matrix(self, phase_difference):
        phase_difference = np.deg2rad(phase_difference)
        retarder_matrix = np.asarray(
            [
                [1, 0, 0, 0],
                [
                    0,
                    cos(2 * self._THETA_POL) ** 2
                    + sin(2 * self._THETA_POL) ** 2 * cos(phase_difference),
                    cos(2 * self._THETA_POL)
                    * sin(2 * self._THETA_POL)
                    * (1 - cos(phase_difference)),
                    -sin(2 * self._THETA_POL) * sin(phase_difference),
                ],
                [
                    0,
                    cos(2 * self._THETA_POL)
                    + sin(2 * self._THETA_POL) * (1 - cos(phase_difference)),
                    sin(2 * self._THETA_POL) ** 2
                    + cos(2 * self._THETA_POL) ** 2 * cos(phase_difference),
                    cos(2 * self._THETA_POL) * sin(phase_difference),
                ],
                [
                    0,
                    sin(2 * self._THETA_POL) * sin(phase_difference),
                    -cos(2 * self._THETA_POL) * sin(phase_difference),
                    cos(phase_difference),
                ],
            ]
        )

        return retarder_matrix


class Concrete_SPARC4_Spectral_Response_1(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 1."""

    _CHANNEL_ID = 1


class Concrete_SPARC4_Spectral_Response_2(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 2."""

    _CHANNEL_ID = 2


class Concrete_SPARC4_Spectral_Response_3(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 3."""

    _CHANNEL_ID = 3


class Concrete_SPARC4_Spectral_Response_4(Abstract_SPARC4_Spectral_Response):

    """Concrete SPARC4 spectral response of the channel 4."""

    _CHANNEL_ID = 4
