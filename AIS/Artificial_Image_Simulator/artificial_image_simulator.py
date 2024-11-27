import datetime
import os
from math import sqrt
from random import randint, uniform

import astropy.io.fits as fits
import numpy as np
from numpy import ndarray

from AIS.Background_Image import Background_Image
from AIS.Header import Header
from AIS.Noise import Noise
from AIS.Point_Spread_Function import Point_Spread_Function
from AIS.Spectral_Energy_Distribution import Sky, Source
from AIS.Spectral_Response import Atmosphere, Channel, Telescope

__all__ = ["Artificial_Image_Simulator"]


class Artificial_Image_Simulator:

    def __init__(
        self,
        ccd_operation_mode: dict[str, int | float | str],
        channel_id: int,
        ccd_temperature: float | int,
    ) -> None:
        """Initialize the class.

        Parameters
        ----------
        ccd_operation_mode : dict[str,int | float | str]
            A python dictionary of the CCD operation mode.
            The allowed keywords are: em_mode, em_gain, preamp readout,
            binn, t_exp, and image_size.
        channel_id : [1, 2, 3, 4]
            Channel id.
        ccd_temperature : float | int
            CCD temperature.

        Example
        -------

        ```
        ccd_operation_mode = {
            'em_mode': 'Conv',
            'em_gain': 1,
            'preamp': 1,
            'readout': 1,
            'binn': 1,
            't_exp': 1,
            'image_size': 100
        }
        ais = Artificial_Image_Simulator(ccd_operation_mode,
                                        channel_id=1, ccd_temperature=-70)
        ```
        """

        self.ccd_operation_mode = ccd_operation_mode
        self.channel_id = channel_id
        self.ccd_temperature = ccd_temperature
        self._verify_ccd_operation_mode()
        self.BGI_obj = Background_Image(
            self.ccd_operation_mode, channel_id, ccd_temperature
        )
        self.SRC_obj = Source()  # TODO isto eh injeção de dependencia
        self.SKY_obj = Sky()
        self.ATM_obj = Atmosphere()
        self.TEL_obj = Telescope()
        self.CHNNL_obj = Channel(channel_id)
        self.PSF_obj = Point_Spread_Function(self.ccd_operation_mode, channel_id)
        self.HDR_obj = Header(self.ccd_operation_mode, ccd_temperature, channel_id)
        return

    @staticmethod
    def _verify_var_in_interval(var, var_name, var_min=0, var_max=2**32) -> None:
        if not var_min <= var <= var_max:
            raise ValueError(
                f"The {var_name} must be in the interval [{var_min},{var_max}]: {var}."
            )

    def _verify_ccd_operation_mode(self) -> None:
        dic_keywords_list = [
            "binn",
            "em_gain",
            "em_mode",
            "image_size",
            "preamp",
            "readout",
            "t_exp",
        ]
        for key in self.ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(f"The provided name is not a CCD parameter: {key}")

        keyvalues = list(self.ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError("There is a missing parameter for the CCD operation mode")

        em_mode = self.ccd_operation_mode["em_mode"]
        em_gain = self.ccd_operation_mode["em_gain"]
        readout = self.ccd_operation_mode["readout"]
        preamp = self.ccd_operation_mode["preamp"]
        binn = self.ccd_operation_mode["binn"]
        t_exp = self.ccd_operation_mode["t_exp"]
        image_size = self.ccd_operation_mode["image_size"]

        self._check_var_in_a_list(em_mode, "EM mode", ["Conv", "EM"])
        if em_mode == "Conv":
            if em_gain != 1:
                raise ValueError(
                    f"The EM Gain must be 1 for the Conventional Mode: {em_gain}"
                )
        else:
            self._verify_var_in_interval(em_gain, "EM Gain", 2, 300)
        self._check_var_in_a_list(readout, "readout rate", [0.1, 1, 10, 20, 30])
        self._check_var_in_a_list(preamp, "pre-amplification", [1, 2])
        self._check_var_in_a_list(binn, "binning", [1, 2])
        self._verify_var_in_interval(image_size, "image size", 1, 1024)
        self.image_size = image_size
        self._verify_var_in_interval(t_exp, "exposure time", 1e-5, 84600)
        self.t_exp = t_exp
        self._verify_var_in_interval(self.ccd_temperature, "ccd temperature", -70, -30)
        return

    @staticmethod
    def _check_var_in_a_list(var, var_name, _list) -> None:
        if var not in _list:
            raise ValueError(f"The allowed values for the {var_name} are: {_list}")

    def create_source_sed_blackbody(
        self,
        magnitude: int | float,
        wavelength_interval: tuple = (),
        temperature: int | float = 0,
    ) -> None:
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
        self.wavelength, self.source_sed = self.SRC_obj.calculate_sed_blackbody(
            magnitude,
            wavelength_interval,
            temperature,
        )

    def create_source_sed_spectral_library(
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
        self.wavelength, self.source_sed = self.SRC_obj.calculate_sed_spectral_library(
            magnitude,
            wavelength_interval,
            spectral_type,
        )

    def apply_linear_polarization(
        self, percent_pol: float = 100, pol_angle: float = 0
    ) -> None:
        """Apply a linear polarization to the SED of the object and the sky.

        Parameters
        ----------
        percent_pol: float, optional
            Percentage of polarization.
        pol_angle: float, optional
            Polarization angle in degrees. If the selected polarization mode were
            linear, the polarization angle must be provided.
        """
        self.source_sed = self.SRC_obj.apply_linear_polarization(percent_pol, pol_angle)

    def apply_circular_polarization(
        self, percent_pol: float = 100, orientation: str = "left"
    ) -> None:
        """Apply ciruclar polarization in SED.

        Parameters
        ----------
        percent_pol: float, optional
            Percentage of polarization.
        orientation: ['right', 'left'], optional
            Orientation of the polarization.
        """
        self.source_sed = self.SRC_obj.apply_circular_polarization(
            percent_pol, orientation
        )

    def apply_polarization(self, stokes: list = []) -> None:
        """Apply a generic polarization to the SED.

        Parameters
        ----------
        stokes: list, optional
            A list with the q, u, and v Stokes parameters.
        """
        self.source_sed = self.SRC_obj.apply_polarization(stokes)

    def apply_Serkowski_curve(self, pol_BVRI: dict) -> None:
        """Apply the Serkowski curve to the SED of the star.

        Parameters
        ----------
        pol_BVRI : dict
            A python dictionary containing the polarization values of the filters BVRI
            in percentage.

        Example
        -------
        ```
        pol_BVRI = {'B':8, 'V': 8.5, 'R': 7.2, 'I':6}
        AIS.apply_Serkowski_curve(pol.BVRI)
        ```
        """
        self.source_sed = self.SRC_obj.apply_Serkowski_curve(pol_BVRI)
        return

    def write_source_sed(self, wavelength: ndarray, sed: ndarray) -> None:
        """Write a customized source SED into the class.

        Parameters
        ----------
        wavelength: array like
            The wavelength interval of the source in nm.
        sed: array like
            Spectral Energy Distribution of the source in W/(m2.m).
        """
        n = sed.shape
        self.source_sed = sed
        if len(n) == 1:
            self.source_sed = np.zeros((4, n[0]), dtype=np.float64)
            self.source_sed[0] = sed
        self.wavelength = wavelength
        self.SRC_obj.write_source_sed(wavelength, sed)
        return

    def print_available_spectral_types(self) -> None:
        """Print the available spectral types."""
        self.SRC_obj.print_available_spectral_types()
        return

    def create_sky_sed(self, moon_phase: str) -> None:
        """Create the Spectral Energy Distribution of the sky.

        Parameters
        ----------
        moon_phase : ['new', 'first quarter', 'third quarter', 'full']
            The phase of the moon.
        """
        self.sky_sed = self.SKY_obj.calculate_sed(moon_phase, self.wavelength)

    def apply_atmosphere_spectral_response(
        self, air_mass: int | float = 1.0, sky_condition: str = "photometric"
    ) -> None:
        """Apply the atmosphere spectral response.

        This functions applies the atmosphere spectral response on the
        Spectral Energy Distribution of the source and the sky.

        Parameters
        ----------
        air_mass: 1.0, optional
            The air mass in the light path.
        sky_condition: ["photometric", "regular", "good"], optional
            The sky condition. According to the value provided for this variable,
            a different extinction coeficient for the atmosphere will be selected.
        """

        self.source_sed = self.ATM_obj.apply_spectral_response(
            self.wavelength, self.source_sed, air_mass, sky_condition
        )

    def apply_telescope_spectral_response(self, date="20230329") -> None:
        """Apply the telescope spectral response.

        This functions applies the telescope spectral response to the
        Spectral Energy Distribution of the source and the sky.

        Parameters
        ----------
        date: ['20230329', '20230328'], optional
            Date of the transmission curve of the telescope.
        """
        self.source_sed = self.TEL_obj.apply_spectral_response(
            self.wavelength, self.source_sed, date
        )
        self.sky_sed = self.TEL_obj.apply_spectral_response(
            self.wavelength, self.sky_sed, date
        )

    def apply_sparc4_spectral_response(
        self,
        acquisition_mode: str,
        calibration_wheel: str = "",
        retarder_waveplate: str = "half",
        retarder_waveplate_angle: int | float = 0,
    ) -> None:
        """Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response to the
        Spectral Energy Distribution of the source and the sky.

        Parameters
        ----------
        acquisition_mode: ["photometry", "polarimetry"]
            The acquisition mode of the sparc4.

        calibration_wheel: ["polarizer", "ideal-polarizer", "ideal-depolarizer"], optional.
            The optical component of the calibration wheel.
            This parameter provides to the user the options of using the real or the
            ideal versions of the optical components of the calibration wheel.
            Based on the provided value, AIS will apply the correspondent Stoke matrix.
            It should be highlighted that the transmission of the optical component
            still be applied.

        retarder_waveplate: ["half", "quarter"], optional
            The waveplate for polarimetric measurements.
            This parameter is used only if the acquisition_mode is 'polarimetry'.

        retarder_waveplate_angle: int | float, optional
            The angle of the retarder waveplate in degrees.
            If the acquisition mode of SPARC4 is polarimetry,
            the retarder waveplate angle should be provided.
        """
        self.CHNNL_obj.write_sparc4_operation_mode(
            acquisition_mode,
            calibration_wheel,
            retarder_waveplate,
            retarder_waveplate_angle,
        )
        self.source_sed = self.CHNNL_obj.apply_spectral_response(
            self.source_sed, self.wavelength
        )
        self.sky_sed = self.CHNNL_obj.apply_spectral_response(
            self.sky_sed, self.wavelength
        )

    def _integrate_sed(self) -> None:
        self.wavelength *= 1e-9
        self.sky_photons_per_second = np.clip(
            np.trapezoid(self.sky_sed, self.wavelength), 0, None
        )
        self.star_photons_per_second = np.clip(
            np.trapezoid(self.source_sed, self.wavelength), 0, None
        )

    def calculate_exposure_time(self, snr: float = 100, seeing: float = 1.5) -> float:
        """Calculate the exposure time required for a given Signal-to-Noise Ratio.

        Parameters
        ----------
        snr: float, optional
            Signal-to-noise ratio. Defaults to 100.
        seeing: float, optional
            Seeing of the object. Defaults to 1.5.

        Returns
        -------
            float
                The minimum exposure time required to meet the provided SNR.
        """
        n_pix = self.PSF_obj.calculate_npix_star(seeing)
        noise_obj = Noise(self.channel_id)
        dark_noise = noise_obj.calculate_dark_current(self.ccd_temperature)
        read_noise = noise_obj.calculate_read_noise(self.ccd_operation_mode)

        noise_factor = 1.0
        em_gain = 1.0
        binn = self.ccd_operation_mode["binn"]
        if self.ccd_operation_mode["em_mode"] == "EM":
            noise_factor = 1.4
            em_gain = self.ccd_operation_mode["em_gain"]

        self._integrate_sed()
        a = self.star_photons_per_second**2
        b = (
            snr**2
            * noise_factor**2
            * (
                self.star_photons_per_second
                + n_pix * (self.sky_photons_per_second + dark_noise)
            )
        )
        c = snr**2 * n_pix * (read_noise / em_gain / binn) ** 2

        texp_list = self._solve_second_degree_eq(a, -b, -c)
        min_t_exp = min([texp for texp in texp_list if texp > 0])

        return min_t_exp

    def create_artificial_image(
        self,
        image_path: str,
        star_coordinates: tuple,
        seeing: float = 1,
        seed: float = 1,
    ):
        """Create an artificial image.

        This function creates a FITS file of an artificial image, similator to those
        acquired using the SPARC4 scientific cameras.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        star_coordinates : tuple
            The coordinates in pixels of the star in the image.
        seeing : float, optional
            The seeing of the star.
        seed: float, optinal.
            The seed used to create the image. Default to 1.

        """
        self._integrate_sed()
        ord_ray, extra_ord_ray = self.star_photons_per_second, 0
        if type(self.star_photons_per_second) != np.float64:
            self.sky_photons_per_second = sum(self.sky_photons_per_second)
            ord_ray, extra_ord_ray = (
                self.star_photons_per_second[0],
                self.star_photons_per_second[1],
            )
        background = self.BGI_obj.create_sky_background(
            self.sky_photons_per_second, seed
        )
        star_psf = self.PSF_obj.create_star_image(
            star_coordinates, ord_ray, extra_ord_ray, seeing, seed
        )
        self._create_image_name(image_path)
        header = self.HDR_obj.create_header()
        header["OBSTYPE"] = "OBJECT"
        header["FILENAME"] = self.image_name
        header["SHUTTER"] = "OPEN"

        image = background + star_psf
        file = os.path.join(image_path, self.image_name)
        fits.writeto(
            file,
            image,
            header=header,
        )

    def create_background_image(
        self, image_path: str, images: int = 1, seed: float = 1
    ) -> None:
        """Create a background image.

        This function creates an image of the sky background.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        images : int, optional
            The number of images to be created. The default is 1.
        seed: float, optinal.
            The seed used to create the image. Default to 1.
        """

        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            self._integrate_sed()
            background = self.BGI_obj.create_sky_background(
                self.sky_photons_per_second, seed
            )
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header["FILENAME"] = self.image_name
            header["OBSTYPE"] = "FLAT"  # ?
            header["SHUTTER"] = "OPEN"

            file = os.path.join(image_path, self.image_name)

            fits.writeto(file, background, header=header)

    def create_dark_image(
        self, image_path: str, images: int = 1, seed: float = 1
    ) -> None:
        """Create a dark image.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        images : int, optional
            The number of dark images to be created. The default is 1.
        seed: float, optinal.
            The seed used to create the image. Default to 1.
        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            dark_image = self.BGI_obj.create_dark_background(seed)
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header["FILENAME"] = self.image_name
            header["OBSTYPE"] = "DARK"
            image_name = os.path.join(image_path, self.image_name)

            fits.writeto(
                image_name,
                dark_image,
                header=header,
            )

    def create_bias_image(
        self, image_path: str, images: int = 1, seed: float = 1
    ) -> None:
        """Create a bias image.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        images : int, optional
            The number of bias images to be created. The default is 1.
        seed: float, optinal.
            The seed used to create the image. Default to 1.
        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            bias = self.BGI_obj.create_bias_background(seed)
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header["EXPTIME"] = 1e-5
            header["FILENAME"] = self.image_name
            header["OBSTYPE"] = "ZERO"

            image_name = os.path.join(image_path, self.image_name)

            fits.writeto(
                image_name,
                bias,
                header=header,
            )

    def create_random_image(
        self, image_path, number_stars=10, seeing: float = 1, seed: float = 1
    ) -> None:
        """Create a random star image.

        This function creates an artificial image with a set of random stars.
        The number of star created by the function can be provided to the class.
        Otherwise, the number 10 will be used.
        Then, all the star images will be sumed with the background image.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        number_stars: int, optional
            The number of stars in the image.
        seeing : float, optional
            The seeing of the star.
        seed: float, optinal.
            The seed used to create the image. Default to 1.
        """
        self._integrate_sed()
        ord_ray, extra_ord_ray = self.star_photons_per_second, 0
        if type(self.star_photons_per_second) != np.float64:
            self.sky_photons_per_second = sum(self.sky_photons_per_second)
            ord_ray, extra_ord_ray = self.star_photons_per_second

        image_size = self.ccd_operation_mode["image_size"]
        star_image = np.zeros((image_size, image_size))
        for _ in range(number_stars):
            x_coord = randint(0, image_size)
            y_coord = randint(0, image_size)
            ordinary_ray = uniform(0, ord_ray)
            extra_ordinary_ray = uniform(0, extra_ord_ray)
            star_image += self.PSF_obj.create_star_image(
                (x_coord, y_coord),
                ordinary_ray,
                extra_ordinary_ray,
                seeing,
            )

        background = self.BGI_obj.create_sky_background(self.sky_photons_per_second)
        self._create_image_name(image_path)
        header = self.HDR_obj.create_header()
        header["OBSTYPE"] = "OBJECT"
        header["FILENAME"] = self.image_name
        header["SHUTTER"] = "OPEN"
        file = os.path.join(image_path, self.image_name)
        fits.writeto(
            file,
            background + star_image,
            header=header,
        )

    def create_flat_image(
        self, image_path: str, images: int = 1, seed: float = 1
    ) -> None:
        """Create a flat image.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.
        images : int, optional
            The number of flat images to be created. The default is 1.
        seed: float, optinal.
            The seed used to create the image. Default to 1.
        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            flat_image = self.BGI_obj.create_flat_background(seed)
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header["FILENAME"] = self.image_name
            header["OBSTYPE"] = "FLAT"
            header["SHUTTER"] = "OPEN"

            image_name = os.path.join(image_path, self.image_name)

            fits.writeto(
                image_name,
                flat_image,
                header=header,
            )

    def _find_image_index(self, image_path) -> int:
        index = 0
        for file in os.listdir(image_path):
            if file.endswith(".fits"):
                new_index = int(file.split("_")[-1][:-5])
                if new_index > index:
                    index = new_index
        return index

    def _create_image_name(self, image_path) -> None:
        now = datetime.datetime.now()
        index = self._find_image_index(image_path)
        self.image_name = now.strftime(
            f"%Y%m%d_s4c{self.channel_id}_{index + 1:06}.fits"
        )
        return

    @staticmethod
    def _solve_second_degree_eq(a, b, c) -> list[float]:
        delta = (b**2) - (4 * a * c)
        x1 = (-b - sqrt(delta)) / (2 * a)
        x2 = (-b + sqrt(delta)) / (2 * a)
        return [x1, x2]
