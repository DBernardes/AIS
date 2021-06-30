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
import sys
from random import randint
from sys import exit

import numpy as np
import pandas as pd

sys.path.append("..")
import astropy.io.fits as fits
from Atmosphere_Spectral_Response import Atmosphere_Spectral_Response
from Background_Image import Background_Image
from Channel_Creator import (
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)
from Header import Header
from Point_Spread_Function import Point_Spread_Function
from Spectrum_Calculation import Spectrum_Calculation
from Telescope_Spectral_Response import Telescope_Spectral_Response


class Artificial_Image_Simulator:
    """Create an image cube with the star flux distribution.

    Parameters
    ----------
    ccd_operation_mode: dictionary
        A python dictionary with the CCD operation mode. The allowed keywords
        values for the dictionary are

        * em_mode: {0, 1}

           Use the 0 for the Conventional Mode and 1 for the EM Mode

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

    gaussian_stddev: int, optional
        Number of pixels of the gaussian standard deviation

    star_coordinates: tuple, optional
        XY star coordinates in the image

    bias_level: int, optional
        Bias level, in ADU, of the image

    image_dir: str, optional
        Directory where the image should be saved

    star_wavelength_interval: (350, 1100, 50), optional
        Wavelength interval of the star for the calculation of the specific flux.

    star_temperature: 5700, optional
        Temperature of the star in Kelvin.

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
        gaussian_std=3,
        star_coordinates=(100, 100),
        bias_level=500,
        sparc4_operation_mode="phot",
        image_dir="",
        star_wavelength_interval=(350, 1150, 50),
        star_temperature=5700,
    ):
        """Initialize the class."""

        if channel in [1, 2, 3, 4]:
            self.channel = channel
        else:
            raise ValueError(
                "There is no camera with the provided" + f" channel: {channel}"
            )

        if type(gaussian_std) is not int:
            raise ValueError(
                "The gaussian standard deviation must be"
                + f"an integer: {gaussian_std}"
            )
        elif gaussian_std <= 0:
            raise ValueError(
                r"The gaussian standard deviation must be greater"
                + f"than zero: {gaussian_std}"
            )
        else:
            self.gaussian_std = gaussian_std

        for coord in star_coordinates:
            if type(coord) is not int:
                raise ValueError(f"The star coordinates must be an integer: {coord}")
            elif coord <= 0:
                raise ValueError(
                    f"The star coordinates must be greater than zero: {coord}"
                )
            else:
                self.star_coordinates = star_coordinates

        if type(bias_level) is not int:
            raise ValueError(f"The bias level must be an integer: {bias_level}")
        elif bias_level <= 0:
            raise ValueError(f"The bias level must be positive: {bias_level}")
        else:
            self.bias_level = bias_level

        if type(sparc4_operation_mode) is not str:
            raise ValueError(
                r"The SPARC4 operation mode must be a string: "
                + f"{sparc4_operation_mode}"
            )
        elif sparc4_operation_mode not in ["phot", "pol"]:
            raise ValueError(
                'The SPARC4 operation mode must be "phot" '
                + f'or "pol": {sparc4_operation_mode}'
            )
        else:
            self.sparc4_operation_mode = sparc4_operation_mode

        if type(image_dir) is not str:
            raise ValueError(f"The directory path must be a string: {image_dir}")
        else:
            self.image_dir = image_dir

        for wavelength in star_wavelength_interval:
            if type(wavelength) is not int:
                raise ValueError(
                    f"The star waveelength interval must be an integer: {wavelength}"
                )
            elif wavelength <= 0:
                raise ValueError(
                    f"The star waveelength interval must be positive: {wavelength}"
                )
            else:
                self.star_wavelength_interval = star_wavelength_interval

        if type(star_temperature) not in [int, float]:
            raise ValueError(
                "The star temperature must be" + f"a number: {star_temperature}"
            )
        elif star_temperature <= 0:
            raise ValueError(
                r"The star temperature must be greater"
                + f"than zero: {star_temperature}"
            )
        else:
            self.star_temperature = star_temperature

        self.ccd_operation_mode = ccd_operation_mode
        self._verify_ccd_operation_mode()
        self._configure_gain()
        self._configure_image_name()

        CHC = 0
        if channel == 1:
            CHC = Concrete_Channel_1(sparc4_operation_mode)
        elif channel == 2:
            CHC = Concrete_Channel_2(sparc4_operation_mode)
        elif channel == 3:
            CHC = Concrete_Channel_3(sparc4_operation_mode)
        elif channel == 4:
            CHC = Concrete_Channel_4(sparc4_operation_mode)
        self.CHC = CHC
        self._calculate_dark_current()
        self._calculate_read_noise(ccd_operation_mode)
        self.PSF = Point_Spread_Function(CHC, ccd_operation_mode, self.ccd_gain)
        self.BGI = Background_Image(
            ccd_operation_mode,
            self.ccd_gain,
            self.dark_current,
            self.read_noise,
            self.bias_level,
        )
        self.HDR = Header(ccd_operation_mode, self.ccd_gain, CHC.get_serial_number())

        self.l_init, self.l_final, self.l_step = (
            self.star_wavelength_interval[0],
            self.star_wavelength_interval[1],
            self.star_wavelength_interval[2],
        )
        self.SC = Spectrum_Calculation(
            temperature=self.star_temperature,
            l_init=self.l_init,
            l_final=self.l_final,
            l_step=self.l_step,
        )
        self.TSR = Telescope_Spectral_Response()
        self.ASR = Atmosphere_Spectral_Response()
        self._calculate_sky_specific_flux()
        self._calculate_star_specific_flux()

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
        ccd_operation_mode = self.ccd_operation_mode
        for key in ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(f"The name provided is not a CCD parameter: {key}")

        keyvalues = list(ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError("There is a missing parameter of the CCD operation mode")

        em_mode = ccd_operation_mode["em_mode"]
        em_gain = ccd_operation_mode["em_gain"]
        hss = ccd_operation_mode["hss"]
        preamp = ccd_operation_mode["preamp"]
        binn = ccd_operation_mode["binn"]
        t_exp = ccd_operation_mode["t_exp"]
        ccd_temp = ccd_operation_mode["ccd_temp"]
        image_size = ccd_operation_mode["image_size"]

        if em_mode not in [0, 1]:
            raise ValueError(f"Invalid value for the EM mode: {em_mode}")
        if em_mode == 0:
            if em_gain != 1:
                raise ValueError(
                    "The EM Gain must be 1 for the Conventional" + f" Mode: {em_gain}"
                )
        else:
            if type(em_gain) not in [float, int]:
                raise ValueError(f"The EM gain must be a number: {em_gain}")
            elif em_gain < 2 or em_gain > 300:
                raise ValueError(f"EM gain out of range [2, 300]: {em_gain}")

        if preamp not in [1, 2]:
            raise ValueError(f"Invalid value for the pre-amplification: {preamp}")

        if hss not in [0.1, 1, 10, 20, 30]:
            raise ValueError(f"Invalid value for the Readout rate: {hss}")

        if binn not in [1, 2]:
            raise ValueError(f"Invalid value for the binning: {bin}")

        if type(image_size) is not int:
            raise ValueError(f"The image size must be an integer: {image_size}")
        elif image_size <= 0:
            raise ValueError(f"The image size must be greater than zero: {image_size}")
        else:
            self.image_size = image_size

        if type(t_exp) not in [float, int]:
            raise ValueError(f"The exposure time must be a number: {t_exp}")
        elif ccd_operation_mode["t_exp"] < 1e-5:
            raise ValueError(f"Invalid value for the exposure time: {t_exp}")

        if type(ccd_temp) not in [float, int]:
            raise ValueError(f"The CCD temperature must be a number: {ccd_temp}")
        if ccd_temp < -80 or ccd_temp > 20:
            raise ValueError(f"CCD temperature out of range [-80, 20]: {ccd_temp}")

    def get_channel_ID(self):
        """Return the ID for the respective SPARC4 channel."""
        return self.CHC.get_channel_ID()

    def _calculate_dark_current(self):
        self.dark_current = self.CHC.calculate_dark_current(
            self.ccd_operation_mode["ccd_temp"]
        )

    def _calculate_read_noise(self, ccd_operation_mode):
        self.read_noise = self.CHC.calculate_read_noise(ccd_operation_mode)

    def _calculate_star_specific_flux(self):
        self.star_specific_flux = self.SC.calculate_star_specific_flux()

    def _calculate_sky_specific_flux(self):
        self.sky_specific_flux = self.SC.calculate_sky_specific_flux()

    def apply_atmosphere_spectral_response(self):
        """Apply the atmosphere spectral response.

        This functions applies the atmosphere spectral response on the
        calculated star specific flux.

        """
        star_specific_flux = self.star_specific_flux
        self.star_specific_flux = self.ASR.apply_atmosphere_spectral_response(
            star_specific_flux,
            l_init=self.l_init,
            l_final=self.l_final,
            l_step=self.l_step,
        )

    def apply_telescope_spectral_response(self):
        """Apply the telescope spectral response.

        This functions applies the telescope spectral response on the
        calculated star specific flux.

        """
        star_specific_flux = self.star_specific_flux
        self.star_specific_flux = self.TSR.apply_telescope_spectral_response(
            star_specific_flux,
            l_init=self.l_init,
            l_final=self.l_final,
            l_step=self.l_step,
        )

    def apply_sparc4_spectral_response(self):
        """Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response on the
        calculated star specific flux.
        """

        (
            self.specific_ordinary_ray,
            self.specific_extra_ordinary_ray,
        ) = self.CHC.apply_sparc4_spectral_response(
            self.star_specific_flux,
            l_init=self.l_init,
            l_final=self.l_final,
            l_step=self.l_step,
        )

    def _integrate_specific_fluxes(self):
        """Integrate the star and the sky specific fluxes."""
        try:
            self.ordinary_ray = np.sum(self.specific_ordinary_ray)
            self.extra_ordinary_ray = np.sum(self.specific_extra_ordinary_ray)
        except Exception:
            self.ordinary_ray = np.sum(self.star_specific_flux)
            self.extra_ordinary_ray = 0
        self.sky_flux = np.sum(self.sky_specific_flux)

    def _configure_image_name(self):
        """Create the image name.

        The image name will be created based on the provided information

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        dic = self.ccd_operation_mode
        em_gain = "_G" + str(dic["em_gain"])
        em_mode = "CONV"
        if dic["em_mode"] == 1:
            em_mode = "EM"
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
            if em_mode == 1:
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
            "Read_Noise_Calculation",
            "spreadsheet",
            f"Channel {self.channel}",
            "Read_noise_and_gain_values.csv",
        )
        spreadsheet = pd.read_csv(file_name)
        # print(spreadsheet["Gain"], tab_index), exit()
        self.ccd_gain = float(spreadsheet["Gain"][tab_index])

    def create_artificial_image(self):
        """Create the artificial star image.

        This function will sum the background image with the star SPF image
        to create an artificil image, similar to those acquired by the
        SPARC4 cameras.


        Returns
        -------
        Star Image:
            A FITS file with the calculated artificial image
        """
        self._integrate_specific_fluxes()
        background = self.BGI.create_background_image(self.sky_flux)
        star_PSF = self.PSF.create_star_PSF(
            self.star_coordinates,
            self.gaussian_std,
            self.ordinary_ray,
            self.extra_ordinary_ray,
        )
        header = self.HDR.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + ".fits")

        fits.writeto(
            image_name,
            background + star_PSF,
            overwrite=True,
            header=header,
        )

    def create_background_image(self):
        """Create the background image.

        This function creates the background image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Background Image:
            A FITS file with the calculated background image
        """
        self._integrate_specific_fluxes()
        background = self.BGI.create_background_image(self.sky_flux)
        header = self.HDR.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + "_BG.fits")

        fits.writeto(
            image_name,
            background,
            overwrite=True,
            header=header,
        )

    def create_dark_image(self):
        """Create a dark image.

        This function creates a dark image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Dark Image:
            A FITS file with the calculated dark image
        """
        dark_image = self.BGI.create_dark_image()
        header = self.HDR.create_header()
        image_name = os.path.join(self.image_dir, self.image_name + "_DARK.fits")

        fits.writeto(
            image_name,
            dark_image,
            overwrite=True,
            header=header,
        )

    def create_bias_image(self):
        """Create a bias image.

        This function creates a bias image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Bias Image:
            A FITS file with the calculated bias image
        """
        bias = self.BGI.create_bias_image()
        header = self.HDR.create_header()
        header["EXPOSURE"] = 1e-5

        image_name = os.path.join(self.image_dir, self.image_name + "_BIAS.fits")

        fits.writeto(
            image_name,
            bias,
            overwrite=True,
            header=header,
        )

    def create_random_image(self, n=10):
        """Create a random star image.

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
        random_image = self.BGI.create_background_image(self.sky_flux)
        for i in range(n):
            x_coord = randint(50, self.image_size - 50)
            y_coord = randint(50, self.image_size - 50)
            ordinary_ray = randint(self.ordinary_ray // 2, self.ordinary_ray // 1)
            extra_ordinary_ray = 0
            if self.sparc4_operation_mode == "pol":
                extra_ordinary_ray = ordinary_ray
            gaussian_std = randint(1, self.gaussian_std)
            random_image += self.PSF.create_star_PSF(
                (x_coord, y_coord), gaussian_std, ordinary_ray, extra_ordinary_ray
            )
        header = self.HDR.create_header()
        image_name = os.path.join(self.image_dir, self.image_name + "_RAND.fits")

        fits.writeto(
            image_name,
            random_image,
            overwrite=True,
            header=header,
        )

    def create_flat_image(self):
        """Create a flat image.

        This function creates a flat image, similar to those
        acquired by the SPARC4 cameras.


        Returns
        -------
        Flat Image:
            A FITS file with the calculated flat image
        """
        flat_image = self.BGI.create_flat_image()
        header = self.HDR.create_header()

        image_name = os.path.join(self.image_dir, self.image_name + "_FLAT.fits")

        fits.writeto(
            image_name,
            flat_image,
            overwrite=True,
            header=header,
        )
