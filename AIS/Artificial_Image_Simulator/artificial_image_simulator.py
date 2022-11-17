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
import datetime
import os
from random import randint
from types import UnionType

import astropy.io.fits as fits
import numpy as np
import pandas as pd
from numpy import ndarray


from ..Background_Image import Background_Image
from ..Header import Header
from ..Point_Spread_Function import Point_Spread_Function
from ..Spectral_Response import Atmosphere, Telescope, Channel
from ..Spectral_Energy_Distribution import Source, Sky


class Artificial_Image_Simulator:
    """Artificial Images Simulator class.

    The Artificial Images Simulator (AIS) class was developed to generate artificial star images,
    similar to those images that would be acquired by using the acquisition system of the SPARC4 instrument.
    To accomplish this, the AIS models as star flux as a 2D gaussian distribution.
    Then, the star flux is added to an image with a background level
    given by counts distribution of an image of the SPARC4 cameras, as a function of its operation mode.

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

    Returns
    ------
        image: ndarray
            An image in the FITS format with the star flux distribution

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
            ccd_operation_mode: dict[str, int | float | str]) -> None:
        """Initialize the class."""
        self.ccd_operation_mode = ccd_operation_mode
        self._verify_ccd_operation_mode()
        return

    @staticmethod
    def _verify_type(var, var_name, _type=float | int):
        if not isinstance(var, _type):
            if type(_type) == UnionType:
                raise ValueError(
                    f"The {var_name} must be an integer or a float: {var}."
                )
            elif _type == int:
                raise ValueError(f"The {var_name} must be an integer: {var}.")
            elif _type == str:
                raise ValueError(f"The {var_name} must be a string: {var}.")
            else:
                pass

    @staticmethod
    def _verify_var_in_interval(var, var_name, var_min=0, var_max=2 ** 32):
        if var <= var_min or var >= var_max:
            raise ValueError(
                f"The {var_name} must be in the interval [{var_min},{var_max}]: {var}."
            )

    def _verify_class_parameters(self) -> None:
        self._verify_type(self.channel, "channel", int)
        self._check_var_in_a_list(self.channel, "channel", [1, 2, 3, 4])

        for coord in self.star_coordinates:
            self._verify_type(coord, "star coordinate")
            self._verify_var_in_interval(coord, "star coordinate")
            if coord > self.ccd_operation_mode["image_size"]:
                raise ValueError(
                    f"The star coordinates must be smaller than the image size: {coord}"
                )

        self._verify_type(self.bias_level, "bias_level", int)
        self._verify_var_in_interval(self.bias_level, "bias_level")

        for wavelength in self.wavelength_interval:
            self._verify_type(wavelength, "wavelength")
            self._verify_var_in_interval(wavelength, "wavelength", 350, 1150)

        self._verify_type(self.seeing, "seeing")
        self._verify_var_in_interval(self.seeing, "seeing")

        self._verify_type(self.star_magnitude, "star magnitude")

        self._verify_type(self.image_dir, "image directory", str)

        self._verify_type(self.moon_condition, "moon condition", str)
        self._check_var_in_a_list(
            self.moon_condition, "moon condition", [
                "new", "waxing", "waning", "full"]
        )
        return

    def _verify_ccd_operation_mode(self) -> None:
        dic_keywords_list = [
            "binn",
            "ccd_temp",
            "em_gain",
            "em_mode",
            "image_size",
            "preamp",
            "readout",
            "t_exp",
        ]
        for key in self.ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(
                    f"The name provided is not a CCD parameter: {key}")

        keyvalues = list(self.ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError(
                "There is a missing parameter of the CCD operation mode")

        em_mode = self.ccd_operation_mode["em_mode"]
        em_gain = self.ccd_operation_mode["em_gain"]
        hss = self.ccd_operation_mode["readout"]
        preamp = self.ccd_operation_mode["preamp"]
        binn = self.ccd_operation_mode["binn"]
        t_exp = self.ccd_operation_mode["t_exp"]
        ccd_temp = self.ccd_operation_mode["ccd_temp"]
        image_size = self.ccd_operation_mode["image_size"]

        self._check_var_in_a_list(em_mode, "EM mode", ["Conv", "EM"])
        if em_mode == "Conv":
            if em_gain != 1:
                raise ValueError(
                    f"The EM Gain must be 1 for the Conventional Mode: {em_gain}"
                )
        else:
            self._verify_type(em_gain, "EM gain")
            self._verify_var_in_interval(em_gain, "EM Gain", 2, 300)

        self._check_var_in_a_list(hss, "readout rate", [0.1, 1, 10, 20, 30])

        self._check_var_in_a_list(preamp, "pre-amplification", [1, 2])

        self._check_var_in_a_list(binn, "binning", [1, 2])

        self._verify_type(image_size, "image size", int)
        self._verify_var_in_interval(image_size, "image size", 0, 1025)
        self.image_size = image_size

        self._verify_type(t_exp, "exposure time")
        self._verify_var_in_interval(t_exp, "exposure time", 1e-5, 84600)
        self.t_exp = t_exp

        self._verify_type(ccd_temp, "CCD temperature")
        self._verify_var_in_interval(ccd_temp, "CCD temperature", -80, 20)
        self.ccd_temp = ccd_temp
        return

    def create_source_sed(self, calculation_method: str,
                          magnitude: int | float,
                          wavelength_interval: tuple = (),
                          temperature: int | float = 0,
                          spectral_type: str = ''):
        """Create the Spectral Energy Distribution of the source.

        Parameters
        ----------

        calculation_method : ['blackbody', 'spectral_library']
            The method used to calculate the SED.

        magnitude : int | float
            The magnitude of the astronomical object in the V band.

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
        """
        src = Source()
        self.wavelength, self.source_sed = src.calculate_sed(calculation_method, magnitude,
                                                             wavelength_interval, temperature, spectral_type)

    def print_available_spectral_types(self) -> None:
        """Print the available spectral types."""
        src = Source()
        src.print_available_spectral_types()
        return

    def create_sky_sed(self, moon_phase: str, object_wavelength: ndarray):
        """Create the Spectral Energy Distribution of the sky.

        Parameters
        ----------
            moon_phase : ['new', 'waning', 'waxing', 'full']
                The phase of the moon.

            object_wavelength : ndarray
                The wavelength of the astronomical object in nm.
        """
        sky = Sky()
        self.sky_sed = sky.calculate_sed(moon_phase, object_wavelength)

    def apply_atmosphere_spectral_response(
        self, air_mass: int | float = 1.0, sky_condition: str = "photometric"
    ):
        """Apply the atmosphere spectral response.

        This functions applies the atmosphere spectral response on the
        Spectral Energy Distribution of the source.

        Parameters
        ----------

        air_mass: 1.0, optional
        The air mass in the light path.

        sky_condition: ["photometric", "regular", "good"], optional
            The sky condition. According to the value provided for this variable,
            a different extinction coeficient for the atmosphere will be selected.
        """

        self._verify_type(air_mass, "air mass")
        self._verify_var_in_interval(air_mass, "air mass", -1e-3)

        self._verify_type(sky_condition, "sky condition", str)
        self._check_var_in_a_list(
            sky_condition, "sky condition", ["photometric", "regular", "good"]
        )

        atm = Atmosphere()
        self.source_sed = atm.apply_spectral_response(
            self.source_sed, self.wavelength, air_mass, sky_condition)

    def apply_telescope_spectral_response(self):
        """Apply the telescope spectral response.

        This functions applies the telescope spectral response on the
        Spectral Energy Distribution of the source and the sky.

        """
        tel = Telescope()
        self.source_sed = tel.apply_spectral_response(
            self.source_sed, self.wavelength)
        self.sky_sed = tel.apply_spectral_response(
            self.sky_sed, self.wavelength)

    def apply_sparc4_spectral_response(
        self,
        channel_id: int | float,
        acquisition_mode: str,
        calibration_wheel: str = "empty",
        retarder_waveplate: str = "half"
    ):
        """Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response on the
        Spectral Energy Distribution of the source and the sky.

        Parameters
        ----------

        channel_id: int | float
            The channel of the SPARC4.

        acquisition_mode: ["photometric", "polarimetric"]
            The acquisition mode of the sparc4.

        calibration_wheel: ["polarizer", "depolarizer", "empty"], optional
            The position of the calibration wheel.
            This parameter is used only if the acquisition_mode is 'polarimetric'.

        retarder: ["half", "quarter"], optional
            The waveplate for polarimetric measurements.
            This parameter is used only if the acquisition_mode is 'polarimetric'.

        """
        channel = Channel(channel_id, acquisition_mode,
                          calibration_wheel, retarder_waveplate)
        self.source_sed = channel.apply_spectral_response(
            self.source_sed, self.wavelength)
        self.sky_sed = channel.apply_spectral_response(
            self.sky_sed, self.wavelength)

    def _integrate_sed(self):
        """Integrate the star and the sky SEDs."""
        self.sky_photons_per_second = np.trapz(self.sky_sed)
        self.star_photons_per_second = np.trapz(self.source_sed)

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
        self._configure_image_name()
        self._integrate_specific_photons_per_second()
        background = self.bgi.create_background_image(
            self.sky_photons_per_second)
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
        self._configure_image_name()
        self._integrate_specific_photons_per_second()
        background = self.bgi.create_background_image(
            self.sky_photons_per_second)
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
        self._configure_image_name()
        dark_image = self.bgi.create_dark_image()
        header = self.hdr.create_header()
        image_name = os.path.join(
            self.image_dir, self.image_name + "_DARK.fits")

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
        self._configure_image_name()
        bias = self.bgi.create_bias_image()
        header = self.hdr.create_header()
        header["EXPOSURE"] = 1e-5

        image_name = os.path.join(
            self.image_dir, self.image_name + "_BIAS.fits")

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
        self._configure_image_name()
        self._integrate_specific_photons_per_second()
        random_image = self.bgi.create_background_image(
            self.sky_photons_per_second)
        for i in range(n):
            image_size = self.image_size
            x_coord = randint(0, image_size)
            y_coord = randint(0, image_size)
            ordinary_ray = randint(0, self.star_ordinary_ray // 1)
            extra_ordinary_ray = 0
            random_image += self.psf.create_star_psf(
                (x_coord, y_coord),
                ordinary_ray,
                extra_ordinary_ray,
            )
        header = self.hdr.create_header()
        image_name = os.path.join(
            self.image_dir, self.image_name + "_RAND.fits")

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
        self._configure_image_name()
        flat_image = self.bgi.create_flat_image()
        header = self.hdr.create_header()

        image_name = os.path.join(
            self.image_dir, self.image_name + "_FLAT.fits")

        fits.writeto(
            image_name,
            flat_image,
            overwrite=True,
            header=header,
        )


def _configure_image_name(self) -> None:
    """Create the image name.

    The image name will be created based on the time that the image is created


    """
    now = datetime.datetime.now()
    self.image_name = (
        f"{now.year}{now.month:0>2}{now.day:0>2}T{now.hour:0>2}{now.minute:0>2}{now.second:0>2}"
        + f"{now.microsecond}"[:2]
    )
    return
