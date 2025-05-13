import os
from copy import copy

# from sbpy.calib import vega_fluxd
from math import cos, pi, sin, sqrt, tan
from sys import exit

import numpy as np
import pandas as pd
from numpy import ndarray
from scipy.constants import c, h, k, pi
from scipy.interpolate import interp1d, splev, splrep
from scipy.optimize import curve_fit

__all__ = ["Source", "Sky"]


class Spectral_Energy_Distribution:

    EFFECT_WAVELENGTH = 545  # nm
    TELESCOPE_EFFECTIVE_AREA = 0.804 * pi * 0.8**2  # m2
    S_0 = 3.631e-2  # W/(m.m2)
    BASE_PATH = os.path.join(os.path.dirname(__file__))

    def __init__(self) -> None:
        """Initialize the class.

        Tips: useful links:
            https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
            https://cass.ucsd.edu/archive/physics/ph162/mags.html
        """
        self.sed = np.linspace(100, 1000, 100)
        return

    def calculate_sed(self) -> ndarray:
        """Calculate the Spectral Energy Distribution.

        This method is an abstract method that must be implemented in the child classes.

        Returns
        -------
        ndarray:
            The Spectral Energy Distribution of the object in photons/s/m.
        """
        return np.linspace(100, 1000, 100, dtype=np.float64)

    @staticmethod
    def _interpolate_spectral_distribution(
        wavelength, spectral_response, obj_wavelength
    ) -> ndarray:
        # spl = interp1d(wavelength, spectral_response,
        #               bounds_error=False, fill_value='extrapolate', kind='cubic')
        spl = splrep(wavelength, spectral_response)
        interpolated_spectral_distribution = splev(obj_wavelength, spl)
        return interpolated_spectral_distribution  # spl(obj_wavelength)

    def _calculate_photons_density(self, magnitude) -> float:
        return (
            self.S_0
            * 10 ** (-magnitude / 2.5)
            * self.TELESCOPE_EFFECTIVE_AREA
            * self.EFFECT_WAVELENGTH
            * 1e-9
            / (h * c)
        )


class Source(Spectral_Energy_Distribution):

    # Kitchin, C. R., 2003. "Astrophysical Techniques"
    effect_wl = {"B": 436e-9, "V": 545e-9, "R": 638e-9, "I": 797e-9}

    def __init__(self) -> None:
        """Initialize the class.

        Example
        -------

        ```
        src = Source()
        wv, sed = src.calculate_sed(calculation_method='blackbody',
                                magnitude=12,
                                wavelength_interval=(400, 1100, 100),
                                temperature=5700)
        ```
        """
        self.SPECTRAL_LIB_PATH = os.path.join(self.BASE_PATH, "Spectral_Library")
        self.pol_BVRI = dict.fromkeys(["B", "V", "R", "I"], None)

        return

    def calculate_sed_blackbody(
        self,
        magnitude: int | float,
        wavelength_interval: tuple = (),
        temperature: int | float = 0,
    ) -> tuple[ndarray]:
        """Calculate the star SED based on the balckbody distribution.


        Parameters
        ----------

        magnitude : int | float
            The magnitude of the astronomical object in the V band.
            The magnitude is used to calculate the effective flux of
            the astronomical object.

        wavelength_interval : tuple, optional
            The wavelength interval, in nm, of the astronomical object.
            This parameter must be a tuple with three elements, where the first
            element is the initial wavelength, the second element is the final
            wavelength and the third element is the number of elements in the array.

        temperature : int | float, optional
            The blackbody temperature of the astronomical object in Kelvin.

        Returns
        -------
            ndarray:
                The wavelength of the astronomical object in nm.
            ndarray:
                The SED of the astronomical object in photons/m/s.
        """

        wavelength = np.linspace(*wavelength_interval, dtype=np.float64)
        sed = self._calculate_sed_blackbody(wavelength, temperature)
        normalization_flux = self._interpolate_spectral_distribution(
            wavelength, sed, self.EFFECT_WAVELENGTH
        )
        sed /= normalization_flux
        effective_flux = self._calculate_photons_density(magnitude)
        self.sed = np.zeros((4, len(wavelength)), dtype=np.float64)
        self.sed[0] = sed * effective_flux
        self.wavelength = wavelength

        return self.wavelength, self.sed

    def calculate_sed_spectral_library(
        self,
        magnitude: int | float,
        wavelength_interval: tuple = (),
        spectral_type: str = "",
    ) -> tuple[ndarray]:
        """Calculate the star SED based on a library of spectral standard stars.

        The spectral response and the wavelength of the object are obtained using a
        library of spectral types. These spectrums are taken from the Library of
        Stellar Spectrum of ESO, and they can be found at:
        https://www.eso.org/sci/facilities/paranal/decommissioned/isaac/tools/lib.html.
        The level of the spectral response is adjusted using the magnitude of the
        object in the V band.

        Parameters
        ----------

        magnitude : int | float
            The magnitude of the astronomical object in the V band.
            The magnitude is used to calculate the effective flux of
            the astronomical object.

        wavelength_interval : tuple, optional
            The wavelength interval, in nm, of the astronomical object.
            This parameter must be a tuple with three elements, where the first
            element is the initial wavelength, the second element is the final
            wavelength and the third element is the number of elements in the array.

        spectral_type : str, optional
            The spectral type of the star that will be used to calculate the SED.
            This parameter is used only if the calculation_method is 'spectral_standard'.
            The available spectral types can be found using
            the `print_available_spectral_types()` method.

        Returns
        -------
            ndarray:
                The wavelength of the astronomical object in nm.
            ndarray:
                The SED of the astronomical object in photons/m/s.
        """

        wavelength = np.linspace(*wavelength_interval, dtype=np.float64)
        lib_wavelength, sed = self._read_spectral_library(spectral_type)
        sed = self._interpolate_spectral_distribution(lib_wavelength, sed, wavelength)

        effective_flux = self._calculate_photons_density(magnitude)
        self.sed = np.zeros((4, len(wavelength)), dtype=np.float64)
        self.sed[0] = sed * effective_flux
        self.wavelength = wavelength

        return self.wavelength, self.sed

    def write_source_sed(self, wavelength: ndarray, sed: ndarray) -> None:
        """Write the wavelength and SED of the source.

        Parameters
        ----------
        wavelength : ndarray
            source wavelength interval in nm.
        sed : ndarray
            Spectral Energy Distribution of the source.
        """

        self.wavelength = wavelength
        self.sed = sed
        return

    def apply_linear_polarization(
        self, percent_pol: float = 100, pol_angle: float = 0
    ) -> ndarray:
        """Apply a linear polarization to the SED.

        Parameters
        ----------
        percent_pol: float, optional
            Percentage of polarization.
        pol_angle: float, optional
            Polarization angle in degrees.
            If the selected polarization mode were linear,
            the polarization angle must be provided.

        Returns
        -------
        sed: ndarray
            Polarized SED.
        """
        if not 0 < percent_pol <= 100:
            raise ValueError(
                f"The percentage of polarization must be in the interval of 0 up to 100: {percent_pol}"
            )
        if not 0 <= pol_angle <= 180:
            raise ValueError(
                f"The polarization angle must be in the interval of 0 up to 180: {pol_angle}"
            )
        theta = np.deg2rad(pol_angle)
        percent_pol /= 100
        self.sed[1] = self.sed[0] * percent_pol * cos(2 * theta)
        self.sed[2] = self.sed[0] * percent_pol * sin(2 * theta)

        # theta = np.deg2rad(pol_angle)
        # percent_pol /= 100
        # tan_value = tan(2 * theta)
        # self.sed[1] = self.sed[0] * percent_pol / sqrt(1 + tan_value**2)
        # self.sed[2] = self.sed[1] * tan_value

        return self.sed

    def apply_circular_polarization(
        self, percent_pol: float = 100, orientation: str = "left"
    ) -> ndarray:
        """Apply a circular polarization to the SED.

        Parameters
        ----------
        percent_pol: float, optional
            Percentage of polarization.
        orientation: ['right', 'left'], optional
            Orientation of the polarization.

        Returns
        -------
        ndarray:
            Polarized SED.
        """
        if not 0 < percent_pol <= 100:
            raise ValueError(
                f"The percentage of polarization must be in the interval of 0 up to 100: {percent_pol}"
            )
        if orientation not in ["right", "left"]:
            raise ValueError(
                f'The value for the orientation must be "left" or "right": {orientation}'
            )

        self.sed[3] = percent_pol * self.sed[0] / 100
        if orientation == "right":
            self.sed[3] *= -1

        return self.sed

    def apply_polarization(self, stokes: list = []) -> ndarray:
        """Apply a generic polarization to the SED.

        Parameters
        ----------
        stokes: list, optional.
            A list of the q, u, and v Stokes parameters.

        Returns
        --------
        ndarray:
            polarized SED, adjusted for the provided Stokes parameters.

        Raises
        -------
        ValueError:
            The variable stokes must have all the three Stokes parameters.
        ValueError:
            The provided values must be in the interval of -1 up to 1.
        ValueError:
            The quadratic sum of the Stokes parameters must be equal or smaller than 1.
        """
        if stokes == []:
            return self.sed
        else:
            stokes = np.asarray(stokes)
            quadratic_sum = sqrt(sum(stokes**2))
            if len(stokes) != 3:
                raise ValueError(
                    f"A wrong value has been provided for the Stokes parameters: {stokes}"
                )
            elif max(stokes) > 1 or min(stokes) < -1:
                raise ValueError(
                    f"The Stokes parameters must be in the interval -1 up to 1: {stokes}"
                )
            elif quadratic_sum > 1:
                raise ValueError(
                    f"The quadratic sum of the Stokes parameters must be equal or smaller than 1: {stokes}"
                )
        I = self.sed[0]
        q, u, v = stokes
        self.sed[1] = I * q
        self.sed[2] = I * u
        self.sed[3] = I * v

        return self.sed

    def apply_Serkowski_curve(self, pol_BVRI: dict, PA: float = 0) -> ndarray:
        """Apply the Serkowski curve to the SED of the star.

        Parameters
        ----------
        pol_BVRI : dict
            A python dictionary containing the polarization values
            of the filters BVRI in percentage.

        PA: float
            Polarization angle in degrees.

        Returns
        -------
        ndarray
            The SED of the star with the q and u Stokes parameters calculated
            according to the Serkowski curve.
        """
        PA = np.deg2rad(PA)
        self._verify_pol_BVRI(pol_BVRI)
        popt = self._adjust_Serkowski_curve(pol_BVRI)
        percent_pol = self._Serkowski_curve(self.wavelength * 1e-9, *popt) / 100

        # self.sed[1] = q_Stokes * self.sed[0]
        # self.sed[2] = q_Stokes * tan(2 * PA)

        self.sed[1] = self.sed[0] * percent_pol * cos(2 * PA)
        self.sed[2] = self.sed[0] * percent_pol * sin(2 * PA)
        return self.sed

    def _adjust_Serkowski_curve(self, pol_BVRI):
        popt, _ = curve_fit(
            self._Serkowski_curve,
            list(self.effect_wl.values()),
            list(pol_BVRI.values()),
            p0=(10, 1.15, 500e-9),
        )
        return popt

    @staticmethod
    def _calculate_sed_blackbody(wavelength, temperature) -> float:
        wavelength = wavelength.copy() * 1e-9
        numerator = 2 * h * c**2 * pi
        denominator = wavelength**5 * (
            np.exp((h * c) / (wavelength * k * temperature)) - 1
        )
        return numerator / denominator

    def _read_spectral_library(self, spectral_type) -> tuple[ndarray]:
        spectral_type = spectral_type.lower()
        path = os.path.join(self.SPECTRAL_LIB_PATH, "uk" + spectral_type + ".csv")
        try:
            ss = pd.read_csv(path, dtype=np.float64)
            return ss["wavelength (nm)"], ss["flux (F_lambda)"]
        except FileNotFoundError:
            print(f"\nThe spectral type {spectral_type} is not available.")
            self.print_available_spectral_types()
            raise FileNotFoundError

    def print_available_spectral_types(self) -> None:
        """Print the available spectral types."""
        spec_types = os.listdir(self.SPECTRAL_LIB_PATH)
        print("\nAvailable spectral types:")
        print("-------------------------\n")
        spec_types = [spec_type.split(".")[0][2:] for spec_type in spec_types]
        print(*spec_types, sep="\n")

    @staticmethod
    def _Serkowski_curve(wavelength, p_max, k, l_max):
        return p_max * np.exp(-k * np.log(l_max / wavelength) ** 2)

    @staticmethod
    def _verify_pol_BVRI(pol_BVRI):
        pol_vals = pol_BVRI.values()
        if None in pol_vals:
            raise ValueError(
                f"The polarization values for all the filters BVRI should be provided: {pol_BVRI}"
            )
        for val in pol_vals:
            if not 0 <= val <= 100:
                raise ValueError(
                    f"The provided polarization values should be in the interval [0, 100]: {val}"
                )

    # def calculate_sed(
    #     self,
    #     calculation_method: str,
    #     magnitude: int | float,
    #     wavelength_interval: tuple = (),
    #     temperature: int | float = 0,
    #     spectral_type: str = "",
    # ) -> tuple[ndarray]:
    #     """Get the Spectral Energy Distribution of the astronomical object.

    #     Parameters
    #     ----------
    #     calculation_method : ['blackbody', 'spectral_library']
    #         The method used to calculate the SED.
    #         If the user wants to use a blackbody SED,
    #         the calculation_method must be 'blackbody'.
    #         If the user wants to use a spectral standard SED,
    #         the calculation_method must be 'spectral_library'.

    #         In the 'blackbody' case, the spectral response of the object is
    #         calculated using the Planck function, given the temperature and the
    #         wavelength interval of the object. In the 'spectral_library' case, the
    #         spectral response and the wavelength of the object are obtained using a
    #         library of spectral types. These spectrums are taken from the Library of
    #         Stellar Spectrum of ESO, and they can be found at:
    #         https://www.eso.org/sci/facilities/paranal/decommissioned/isaac/tools/lib.html.
    #         The level of the spectral response is adjusted using the magnitude of the
    #         object in the V band.

    #     magnitude : int | float
    #         The magnitude of the astronomical object in the V band.
    #         The magnitude is used to calculate the effective flux of
    #         the astronomical object.

    #     wavelength_interval : tuple, optional
    #         The wavelength interval, in nm, of the astronomical object.
    #         This parameter must be a tuple with three elements, where the first
    #         element is the initial wavelength,
    #         the second element is the final wavelength and the third element is the
    #         number of elements in the array.
    #         This parameter is used only if the calculation_method is 'blackbody'.

    #     temperature : int | float, optional
    #         The blackbody temperature of the astronomical object in Kelvin.
    #         This parameter is used only if the calculation_method is 'blackbody'.

    #     spectral_type : str, optional
    #         The spectral type of the star that will be used to calculate the SED.
    #         This parameter is used only if the calculation_method is 'spectral_standard'.
    #         The available spectral types can be found using
    #         the `print_available_spectral_types()` method.

    #     Returns
    #     -------
    #         ndarray:
    #             The wavelength of the astronomical object in nm.
    #         ndarray:
    #             The SED of the astronomical object in photons/m/s.
    #     """

    #     init, final, step = wavelength_interval
    #     wavelength = np.linspace(init, final, step, dtype=np.float64)
    #     if calculation_method == "blackbody":
    #         sed = self._calculate_sed_blackbody(wavelength, temperature)
    #         normalization_flux = self._interpolate_spectral_distribution(
    #             wavelength, sed, self.EFFECT_WAVELENGTH
    #         )
    #         sed /= normalization_flux
    #     elif calculation_method == "spectral_library":
    #         lib_wavelength, sed = self._read_spectral_library(spectral_type.lower())
    #         sed = self._interpolate_spectral_distribution(
    #             lib_wavelength, sed, wavelength
    #         )

    #     else:
    #         raise ValueError(
    #             f"The calculation_method should be 'blackbody' or 'spectral_library': {calculation_method}."
    #         )

    #     n = len(sed)
    #     effective_flux = self._calculate_photons_density(magnitude)
    #     self.sed = np.zeros((4, n), dtype=np.float64)
    #     self.sed[0] = sed * effective_flux
    #     self.wavelength = wavelength

    #     return self.wavelength, self.sed


class Sky(Spectral_Energy_Distribution):

    CSV_FILE = "moon_magnitude.csv"

    def __init__(self) -> None:
        """Initialize the Sky class."""
        return

    def _read_csv(self, file_name, value_name) -> tuple[ndarray]:
        file_name = os.path.join(self.BASE_PATH, file_name)
        ss = pd.read_csv(file_name, dtype=np.float64)
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
        ndarray:
            The wavelength of the sky in nm.
        ndarray:
            The SED of the sky in photons/m/s.
        """
        wavelength, mags = self._read_csv(self.CSV_FILE, moon_phase)
        sed = self._calculate_photons_density(mags)
        temp = self._interpolate_spectral_distribution(
            wavelength, sed, object_wavelength
        )
        sed = np.zeros((4, object_wavelength.shape[0]), dtype=np.float64)
        sed[0] = temp
        return sed

    @staticmethod
    def _interpolate_spectral_distribution(
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
