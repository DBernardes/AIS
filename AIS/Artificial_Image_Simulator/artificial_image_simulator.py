"""
Artificial Images Simulator
============================

The Artificial Images Simulator (AIS) class was developed to generate
artificial star images, similar to those images that would be acquired by
using the acquisition system of the instrument. To accomplish this,
the AIS models as star flux as a 2D gaussian distribution. Then, the star
flux is added to an image with a background level given by counts distribution
of an image of the SPARC4 cameras, as a function of its operation mode.
"""
import os
from random import randint

import astropy.io.fits as fits
import numpy as np
import pandas as pd

from ..Atmosphere_Spectral_Response import Atmosphere_Spectral_Response
from ..Background_Image import Background_Image
from ..Channel_Creator import (
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)
from ..Header import Header
from ..Point_Spread_Function import Point_Spread_Function
from ..Spectrum_Calculation import Spectrum_Calculation
from ..Telescope_Spectral_Response import Telescope_Spectral_Response

# from sys import exit


class Artificial_Image_Simulator:
    """
    Create an image cube with the star flux distribution.

    Parameters
    ----------
    ccd_operation_mode: dictionary
        A python dictionary with the CCD operation mode. The allowed keywords
        values for the dictionary are

        * em_mode: {Conv, EM}

           Use the Conv for the Conventional Mode and EM for the Electron Multiplying mode

        * em_gain: float

           Electron Multiplying gain

        * preamp: {1, 2}

           Pre-amplification

        * hss: {0.1, 1, 10, 20, 30}

           Horizontal Shift Speed (readout rate) in MHz

        * bin: int

           Number of the binned pixels

        * t_exp: float

           Exposure time in seconds

        * ccd_temp: float

            CCD temperature

        * image_size: int, optional

            Image size in pixels

    channel: {1, 2, 3, 4}, optional
        The SPARC4 channel

    star_coordinates: tuple, optional
        XY star coordinates in the image

    bias_level: int, optional
        Bias level, in ADU, of the image

    image_dir: str, optional
        Directory where the image should be saved

    wavelength_interval: (400, 1100, 50), optional
        Wavelength interval of the star for the calculation of the specific flux.

    star_temperature: 5700, optional
        Temperature of the star in Kelvin.

    star_magnitude: 22, optional
        Magnitude of the star.

    seeing: 1.5, optional
        Seeing disc of the observatory.

    sparc4_operation_mode: dictionary
        A python dictionary with the SPARC4 operation mode. The allowed keywords for the
        dictionary are:

        * acquisition_mode: {"photometry", "polarimetry"}

            The acquisition mode of the sparc4.

        * calibration_wheel: {"polarizer", "depolarizer", "empty"}

            The position of the calibration wheel.

        * waveplate: {"half", "quarter"}

            The waveplate for polarimetric measurements.

    air_mass: 1.0, optional
        The air mass in the light path.

    sky_condition: {"photometric", "regular", "good"}
        The condition of the sky at the observaiton night. According to the value provided for this variable,
        a different extinction coeficient for the atmosphere will be selected.

    Yields
    ------
        image cube: array like
            An image cube in the FITS format with the star flux distribution

    Notes
    -----
        Explicar o código; background; passo-a-passo

    Examples
    --------
        Incluir exemplos

    References
    ----------
    .. [#Bernardes_2018] Bernardes, D. V., Martioli, E., and Rodrigues, C. V., “Characterization of the SPARC4 CCDs”, <i>Publications of the Astronomical Society of the Pacific</i>, vol. 130, no. 991, p. 95002, 2018. doi:10.1088/1538-3873/aacb1e.

    """

    def __init__(
        self,
        ccd_operation_mode,
        channel=1,
        star_coordinates=(100, 100),
        bias_level=500,
        sparc4_operation_mode="phot",
        image_dir="",
        wavelength_interval=(400, 1100, 50),
        star_temperature=5700,
        star_magnitude=22,
        seeing=1.5,
        air_mass=1.0,
        sky_condition="good",
    ):
        """Initialize the class."""
        self.ccd_operation_mode = ccd_operation_mode
        self.channel = channel
        self.star_coordinates = star_coordinates
        self.bias_level = bias_level
        self.sparc4_operation_mode = sparc4_operation_mode
        self.image_dir = image_dir
        self.wavelength_interval = wavelength_interval
        self.star_temperature = star_temperature
        self.star_magnitude = star_magnitude
        self.seeing = seeing
        self.air_mass = air_mass
        self.sky_condition = sky_condition

        self._verify_ccd_operation_mode()
        self._verify_class_parameters()

        l_init, l_final, l_step = (
            self.wavelength_interval[0],
            self.wavelength_interval[1],
            self.wavelength_interval[2],
        )

        self.wavelength_interval = range(l_init, l_final + l_step, l_step)
        self.wavelength_len = len(self.wavelength_interval)

        self._configure_gain()
        self._configure_image_name()
        self._initialize_subclasses()
        self._calculate_sky_specific_flux()
        self._calculate_star_specific_flux()

    def _verify_class_parameters(self):
        if self.channel not in [1, 2, 3, 4]:
            raise ValueError(
                "There is no camera with the provided" + f" channel: {self.channel}"
            )

        for coord in self.star_coordinates:
            if not isinstance(coord, int):
                raise ValueError(f"The star coordinates must be an integer: {coord}")
            elif coord <= 0:
                raise ValueError(
                    f"The star coordinates must be greater than zero: {coord}"
                )
            elif coord > self.ccd_operation_mode["image_size"]:
                raise ValueError(
                    f"The star coordinates must be smaller than the image size: {coord}"
                )

        if not isinstance(self.bias_level, int):
            raise ValueError(f"The bias level must be an integer: {self.bias_level}")
        elif self.bias_level <= 0:
            raise ValueError(f"The bias level must be positive: {self.bias_level}")

        if not isinstance(self.sparc4_operation_mode, str):
            raise ValueError(
                r"The SPARC4 operation mode must be a string: "
                + f"{self.sparc4_operation_mode}"
            )
        elif self.sparc4_operation_mode not in ["phot", "pol"]:
            raise ValueError(
                'The SPARC4 operation mode must be "phot" '
                + f'or "pol": {self.sparc4_operation_mode}'
            )

        if not isinstance(self.image_dir, str):
            raise ValueError(f"The directory path must be a string: {self.image_dir}")

        for wavelength in self.wavelength_interval:
            if not isinstance(wavelength, int):
                raise ValueError(
                    f"The wavelength interval must be an integer: {wavelength}"
                )
            elif wavelength <= 0:
                raise ValueError(
                    f"The wavelength interval must be positive: {wavelength}"
                )

        if not isinstance(self.star_temperature, (int, float)):
            raise ValueError(
                "The star temperature must be" + f"a number: {self.star_temperature}"
            )
        elif self.star_temperature <= 0:
            raise ValueError(
                r"The star temperature must be greater"
                + f"than zero: {self.star_temperature}"
            )

        if not isinstance(self.star_magnitude, (int, float)):
            raise ValueError(
                "The star magnitude must be" + f"a number: {self.star_magnitude}"
            )

        if not isinstance(self.seeing, (int, float)):
            raise ValueError("The seeing must be" + f"a number: {self.seeing}")
        elif self.seeing <= 0:
            raise ValueError(
                r"The seeing must be greater" + f"than zero: {self.seeing}"
            )

        if not isinstance(self.air_mass, (int, float)):
            raise ValueError("The air mass must be" + f"a number: {self.air_mass}")
        elif self.air_mass < 0:
            raise ValueError(
                r"The air mass must be greater or equal to" + f"zero: {self.air_mass}"
            )

        if not isinstance(self.sky_condition, str):
            raise ValueError(
                "The sky condition must be" + f"a string: {self.sky_condition}"
            )
        elif self.sky_condition not in ["photometric", "regular", "good"]:
            raise ValueError(
                r"The allowed values for the sky condition are 'photometric', 'regular', or 'good'."
            )

    def _verify_ccd_operation_mode(self):
        """Verify if the provided CCD operation mode is correct."""
        dic_keywords_list = [
            "binn",
            "ccd_temp",
            "em_gain",
            "em_mode",
            "hss",
            "image_size",
            "preamp",
            "t_exp",
        ]
        for key in self.ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(f"The name provided is not a CCD parameter: {key}")

        keyvalues = list(self.ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError("There is a missing parameter of the CCD operation mode")

        em_mode = self.ccd_operation_mode["em_mode"]
        em_gain = self.ccd_operation_mode["em_gain"]
        hss = self.ccd_operation_mode["hss"]
        preamp = self.ccd_operation_mode["preamp"]
        binn = self.ccd_operation_mode["binn"]
        t_exp = self.ccd_operation_mode["t_exp"]
        ccd_temp = self.ccd_operation_mode["ccd_temp"]
        image_size = self.ccd_operation_mode["image_size"]

        if em_mode not in ["Conv", "EM"]:
            raise ValueError(f"Invalid value for the EM mode: {em_mode}")
        if em_mode == "Conv":
            if em_gain != 1:
                raise ValueError(
                    f"The EM Gain must be 1 for the Conventional Mode: {em_gain}"
                )
        else:
            if not isinstance(em_gain, (float, int)):
                raise ValueError(f"The EM gain must be a number: {em_gain}")
            elif em_gain < 2 or em_gain > 300:
                raise ValueError(f"EM gain out of range [2, 300]: {em_gain}")

        if preamp not in [1, 2]:
            raise ValueError(f"Invalid value for the pre-amplification: {preamp}")

        if hss not in [0.1, 1, 10, 20, 30]:
            raise ValueError(f"Invalid value for the Readout rate: {hss}")

        if binn not in [1, 2]:
            raise ValueError(f"Invalid value for the binning: {bin}")

        if not isinstance(image_size, int):
            raise ValueError(f"The image size must be an integer: {image_size}")
        elif image_size <= 0:
            raise ValueError(f"The image size must be greater than zero: {image_size}")
        else:
            self.image_size = image_size

        if not isinstance(t_exp, (float, int)):
            raise ValueError(f"The exposure time must be a number: {t_exp}")
        elif t_exp < 1e-5:
            raise ValueError(f"Invalid value for the exposure time: {t_exp}")

        if not isinstance(ccd_temp, (float, int)):
            raise ValueError(f"The CCD temperature must be a number: {ccd_temp}")
        if ccd_temp < -80 or ccd_temp > 20:
            raise ValueError(f"CCD temperature out of range [-80, 20]: {ccd_temp}")

    def _initialize_subclasses(self):
        self.sc = Spectrum_Calculation(
            star_temperature=self.star_temperature,
            wavelength_interval=self.wavelength_interval,
        )
        self.tsr = Telescope_Spectral_Response()
        self.asr = Atmosphere_Spectral_Response(self.air_mass, self.sky_condition)

        # -------------------------------------------------------------------------------------------

        chc = 0
        if self.channel == 1:
            chc = Concrete_Channel_1(
                self.sparc4_operation_mode, self.wavelength_interval
            )
        elif self.channel == 2:
            chc = Concrete_Channel_2(
                self.sparc4_operation_mode, self.wavelength_interval
            )
        elif self.channel == 3:
            chc = Concrete_Channel_3(
                self.sparc4_operation_mode, self.wavelength_interval
            )
        elif self.channel == 4:
            chc = Concrete_Channel_4(
                self.sparc4_operation_mode, self.wavelength_interval
            )
        self.chc = chc
        self._calculate_dark_current()
        self._calculate_read_noise(self.ccd_operation_mode)
        self.psf = Point_Spread_Function(
            chc, self.ccd_operation_mode, self.ccd_gain, self.seeing
        )
        self.bgi = Background_Image(
            self.ccd_operation_mode,
            self.ccd_gain,
            self.dark_current,
            self.read_noise,
            self.bias_level,
        )
        self.hdr = Header(
            self.ccd_operation_mode, self.ccd_gain, chc.get_serial_number()
        )

    def get_channel_id(self):
        """Return the ID for the respective SPARC4 channel."""
        return self.chc.get_channel_id()

    def _calculate_dark_current(self):
        self.dark_current = self.chc.calculate_dark_current(
            self.ccd_operation_mode["ccd_temp"]
        )

    def _calculate_read_noise(self, ccd_operation_mode):
        self.read_noise = self.chc.calculate_read_noise(ccd_operation_mode)

    def _calculate_star_specific_flux(self):
        self.star_specific_photons_per_second = [
            self.sc.calculate_specific_flux(self.star_magnitude)
        ]

    def _calculate_sky_specific_flux(self):
        self.sky_specific_photons_per_second = [
            self.sc.calculate_specific_flux(self.star_magnitude + 3)
        ]

    def apply_atmosphere_spectral_response(self):
        """
        Apply the atmosphere spectral response.

        This functions applies the atmosphere spectral response on the
        calculated star and sky specific flux.

        """

        self.star_specific_photons_per_second = [
            self.asr.apply_atmosphere_spectral_response(
                self.star_specific_photons_per_second[0],
                wavelength_interval=self.wavelength_interval,
            )
        ]

        self.sky_specific_photons_per_second = [
            self.asr.apply_atmosphere_spectral_response(
                self.sky_specific_photons_per_second[0],
                wavelength_interval=self.wavelength_interval,
            )
        ]

    def apply_telescope_spectral_response(self):
        """
        Apply the telescope spectral response.

        This functions applies the telescope spectral response on the
        calculated star and sky specific flux.

        """

        self.star_specific_photons_per_second = [
            self.tsr.apply_telescope_spectral_response(
                self.star_specific_photons_per_second[0],
                wavelength_interval=self.wavelength_interval,
            )
        ]

        self.sky_specific_photons_per_second = [
            self.tsr.apply_telescope_spectral_response(
                self.sky_specific_photons_per_second[0],
                wavelength_interval=self.wavelength_interval,
            )
        ]

    def apply_sparc4_spectral_response(self):
        """
        Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response on the
        calculated star and sky specific flux.
        """
        self.star_specific_photons_per_second = self.chc.apply_sparc4_spectral_response(
            self.star_specific_photons_per_second[0]
        )
        self.sky_specific_photons_per_second = self.chc.apply_sparc4_spectral_response(
            self.sky_specific_photons_per_second[0]
        )

    def _integrate_specific_fluxes(self):
        """Integrate the star and the sky specific fluxes."""
        self.sky_photons_per_second = 0
        for array in self.sky_specific_photons_per_second:
            self.sky_photons_per_second += np.trapz(array[0, :])
        star_photons_per_second = []
        for array in self.star_specific_photons_per_second:
            star_photons_per_second.append(np.trapz(array[0, :]))
        self.star_ordinary_ray = star_photons_per_second[0]
        self.star_extra_ordinary_ray = 0
        if len(star_photons_per_second) > 1:
            self.star_extra_ordinary_ray = star_photons_per_second[1]

    def _configure_image_name(self):
        """
        Create the image name.

        The image name will be created based on the provided information

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        dic = self.ccd_operation_mode
        em_gain = "_G" + str(dic["em_gain"])
        em_mode = dic["em_mode"].upper()
        hss = "_HSS" + str(dic["hss"])
        preamp = "_PA" + str(dic["preamp"])
        binn = "_B" + str(dic["binn"])
        t_exp = "_TEXP" + str(dic["t_exp"])
        self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain
        if self.sparc4_operation_mode == "pol":
            self.image_name += "_POL"

    def _configure_gain(self):
        """Configure the CCD gain based on its operation mode."""
        ccd_operation_mode = self.ccd_operation_mode
        em_mode = ccd_operation_mode["em_mode"]
        hss = ccd_operation_mode["hss"]
        preamp = ccd_operation_mode["preamp"]
        tab_index = 0
        if hss == 0.1:
            tab_index = 21
        elif hss == 1:
            tab_index = 17
            if em_mode == "EM":
                tab_index = 13
        elif hss == 10:
            tab_index = 9
        elif hss == 20:
            tab_index = 5
        elif hss == 30:
            tab_index = 1
        else:
            raise ValueError(f"Unexpected value for the readout rate: {hss}")
        if preamp == 2:
            tab_index += 2
        file_name = os.path.join(
            "AIS",
            "Read_Noise_Calculation",
            "spreadsheet",
            f"Channel {self.channel}",
            "Read_noise_and_gain_values.csv",
        )
        spreadsheet = pd.read_csv(file_name)
        self.ccd_gain = float(spreadsheet["Gain"][tab_index])

    def create_artificial_image(self):
        """
        Create the artificial star image.

        This function will sum the background image with the star SPF image
        to create an artificil image, similar to those acquired by the
        SPARC4 cameras.


        Returns
        -------
        Star Image:
            A FITS file with the calculated artificial image
        """
        self._integrate_specific_fluxes()
        background = self.bgi.create_background_image(self.sky_photons_per_second)
        star_psf = self.psf.create_star_psf(
            self.star_coordinates,
            self.star_ordinary_ray,
            self.star_extra_ordinary_ray,
        )
        header = self.hdr.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + ".fits")

        fits.writeto(
            image_name,
            background + star_psf,
            overwrite=True,
            header=header,
        )

    def create_background_image(self):
        """
        Create the background image.

        This function creates the background image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Background Image:
            A FITS file with the calculated background image
        """
        self._integrate_specific_fluxes()
        background = self.bgi.create_background_image(self.sky_photons_per_second)
        header = self.hdr.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + "_BG.fits")

        fits.writeto(
            image_name,
            background,
            overwrite=True,
            header=header,
        )

    def create_dark_image(self):
        """
        Create a dark image.

        This function creates a dark image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Dark Image:
            A FITS file with the calculated dark image
        """
        dark_image = self.bgi.create_dark_image()
        header = self.hdr.create_header()
        image_name = os.path.join(self.image_dir, self.image_name + "_DARK.fits")

        fits.writeto(
            image_name,
            dark_image,
            overwrite=True,
            header=header,
        )

    def create_bias_image(self):
        """
        Create a bias image.

        This function creates a bias image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Bias Image:
            A FITS file with the calculated bias image
        """
        bias = self.bgi.create_bias_image()
        header = self.hdr.create_header()
        header["EXPOSURE"] = 1e-5

        image_name = os.path.join(self.image_dir, self.image_name + "_BIAS.fits")

        fits.writeto(
            image_name,
            bias,
            overwrite=True,
            header=header,
        )

    def create_random_image(self, n=10):
        """
        Create a random star image.

        This function creates an artificial image with a set of random stars.
        The number of star created by the function can be provided to the class. Otherwise,
        the number 10 will be used. Then, all the star images will be sumed with the background image.

        Parameters
        ----------

        n: int, optional
            The number of stars in the image.

        Returns
        -------
        Star Image:
            A FITS file with the calculated random artificial image.
        """

        self._integrate_specific_fluxes()
        random_image = self.bgi.create_background_image(self.sky_photons_per_second)
        for i in range(n):
            image_size = self.image_size
            x_coord = randint(0, image_size)
            y_coord = randint(0, image_size)
            ordinary_ray = randint(0, self.star_ordinary_ray // 1)
            extra_ordinary_ray = 0
            if self.sparc4_operation_mode == "pol":
                extra_ordinary_ray = ordinary_ray
            random_image += self.psf.create_star_psf(
                (x_coord, y_coord),
                ordinary_ray,
                extra_ordinary_ray,
            )
        header = self.hdr.create_header()
        image_name = os.path.join(self.image_dir, self.image_name + "_RAND.fits")

        fits.writeto(
            image_name,
            random_image,
            overwrite=True,
            header=header,
        )

    def create_flat_image(self):
        """
        Create a flat image.

        This function creates a flat image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Flat Image:
            A FITS file with the calculated flat image
        """
        flat_image = self.bgi.create_flat_image()
        header = self.hdr.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + "_FLAT.fits")

        fits.writeto(
            image_name,
            flat_image,
            overwrite=True,
            header=header,
        )
