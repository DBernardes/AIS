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
import sys
from random import randint

sys.path.append("..")
import astropy.io.fits as fits
import openpyxl
from ATM_SR import Atmosphere_Spectral_Response
from BGI import Background_Image
from CHC import (
    Concrete_Channel_1,
    Concrete_Channel_2,
    Concrete_Channel_3,
    Concrete_Channel_4,
)
from HDR import Header
from PSF import Point_Spread_Function
from SC import Spectrum_Calculation
from TEL_SR import Telescope_Spectral_Response


class Artificial_Image_Simulator:
    """Create an image cube with the star flux distribution.

    Parameters
    ----------
    star_magitude : float
        Magnitude of the star
    sky_magnitude: float
        Magnitude of the sky
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

    gaussian_stddev: int
        Number of pixels of the gaussian standard deviation

    star_coordinates: tuple
        XY star coordinates in the image

    bias_level: int, optional
        Bias level, in ADU, of the image

    image_dir: str, optional
        Directory where the image should be saved


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
        star_magnitude,
        sky_magnitude,
        ccd_operation_mode,
        channel=1,
        gaussian_std=3,
        star_coordinates=(100, 100),
        bias_level=500,
        sparc4_operation_mode="phot",
        image_dir="",
    ):
        """Initialize the class."""
        if type(star_magnitude) not in [int, float]:
            raise ValueError("The star flux must be a number: " + f"{star_magnitude}")
        elif star_magnitude <= 0:
            raise ValueError(
                f"The star flux must be greater than zero: {star_magnitude}"
            )
        else:
            self.star_magnitude = star_magnitude

        if type(sky_magnitude) not in [int, float]:
            raise ValueError(f"The sky flux must be a number: {sky_magnitude}")
        elif sky_magnitude <= 0:
            raise ValueError(f"The sky flux must be greater than zero: {sky_magnitude}")
        else:
            self.sky_magnitude = sky_magnitude

        if channel in [1, 2, 3, 4]:
            self.channel = channel
        else:
            raise ValueError(
                "There is no camera with the provided" + f"serial number: {channel}"
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
                raise ValueError("The star coordinates must be an integer: {coord}")
            elif coord <= 0:
                raise ValueError(
                    "The star coordinates must be greater than zero: {coord}"
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
            if image_dir != "":
                if "/" not in image_dir[-1]:
                    image_dir += "\\"
            self.image_dir = image_dir

        self._verify_ccd_operation_mode(ccd_operation_mode)
        self._configure_gain(ccd_operation_mode)
        self._configure_image_name(ccd_operation_mode)

        CHC = 0
        ccd_temp = ccd_operation_mode["ccd_temp"]
        if channel == 1:
            CHC = Concrete_Channel_1(
                ccd_temp, sparc4_operation_mode=sparc4_operation_mode
            )
        elif channel == 2:
            CHC = Concrete_Channel_2(
                ccd_temp, sparc4_operation_mode=sparc4_operation_mode
            )
        elif channel == 3:
            CHC = Concrete_Channel_3(
                ccd_temp, sparc4_operation_mode=sparc4_operation_mode
            )
        elif channel == 4:
            CHC = Concrete_Channel_4(
                ccd_temp, sparc4_operation_mode=sparc4_operation_mode
            )
        self.CHC = CHC
        self._calculate_dark_current()
        self._calculate_read_noise(ccd_operation_mode)
        self.PSF = Point_Spread_Function(
            CHC,
            ccd_operation_mode,
            self.ccd_gain,
        )
        self.BGI = Background_Image(
            ccd_operation_mode,
            self.ccd_gain,
            self.dark_current,
            self.read_noise,
            self.bias_level,
        )
        self.HDR = Header(ccd_operation_mode, self.ccd_gain, CHC.get_serial_number())
        self.SC = Spectrum_Calculation()
        self.TSR = Telescope_Spectral_Response()
        self.ASR = Atmosphere_Spectral_Response()
        self._calculate_sky_spectrum()
        self._calculate_star_spectrum()

    def _verify_ccd_operation_mode(self, ccd_operation_mode):
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
            if em_gain not in [float, int]:
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
        self.dark_current = self.CHC.calculate_dark_current()

    def _calculate_read_noise(self, ccd_operation_mode):
        self.read_noise = self.CHC.calculate_read_noise(ccd_operation_mode)

    def _calculate_star_spectrum(self):
        self.star_spectrum = self.SC.calculate_star_spectrum()

    def _calculate_sky_spectrum(self):
        self.sky_spectrum = self.SC.calculate_sky_spectrum()

    def apply_atmosphere_spectral_response(self):
        """Apply the atmosphere spectral response.

        This functions applies the atmosphere spectral response on the
        calculated star spectrum.

        """
        star_spectrum = self.star_spectrum
        self.star_spectrum = self.ASR.apply_atmosphere_spectral_response(star_spectrum)

    def apply_telescope_spectral_response(self):
        """Apply the telescope spectral response.

        This functions applies the telescope spectral response on the
        calculated star spectrum.

        """
        star_spectrum = self.star_spectrum
        self.star_spectrum = self.TSR.apply_telescope_spectral_response(star_spectrum)

    def apply_sparc4_spectral_response(self):
        """Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response on the
        calculated star spectrum.
        """
        star_spectrum = self.star_spectrum
        self.star_spectrum = self.CHC.apply_sparc4_spectral_response(star_spectrum)

    def _integrate_spectruns(self):
        """Integra the star and the sky spectrums."""
        self.star_flux = sum(self.star_spectrum)
        self.sky_flux = sum(self.sky_spectrum)

    def _configure_image_name(self, ccd_operation_mode, include_star_mag=False):
        """Create the image name.

        The image name will be created based on the provided information

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        dic = ccd_operation_mode
        em_gain = "_G" + str(dic["em_gain"])
        em_mode = "CONV"
        if dic["em_mode"] == 1:
            em_mode = "EM"
        hss = "_HSS" + str(dic["hss"])
        preamp = "_PA" + str(dic["preamp"])
        binn = "_B" + str(dic["binn"])
        t_exp = "_TEXP" + str(dic["t_exp"])
        self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain

        if include_star_mag:
            star_flux = "_S" + str(self.star_magnitude)
            self.image_name += star_flux

    def _configure_gain(self, ccd_operation_mode):
        """Configure the CCD gain based on its operation mode."""
        em_mode = ccd_operation_mode["em_mode"]
        hss = ccd_operation_mode["hss"]
        preamp = ccd_operation_mode["preamp"]
        tab_index = 0
        if hss == 0.1:
            tab_index = 23
        elif hss == 1:
            tab_index = 19
            if em_mode == 1:
                tab_index = 15
        elif hss == 10:
            tab_index = 11
        elif hss == 20:
            tab_index = 7
        elif hss == 30:
            tab_index = 3
        else:
            raise ValueError("Unexpected value for the readout rate: {hss}")
        if preamp == 2:
            tab_index += 2

        spreadsheet = openpyxl.load_workbook(
            f"./RNC/spreadsheet/Channel {self.channel}"
            + "/Read_noise_and_gain_values.xlsx"
        ).active
        self.ccd_gain = spreadsheet.cell(tab_index, 5).value

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
        self._integrate_spectruns()
        background = self.BGI.create_background_image(self.sky_flux)
        star_PSF = self.PSF.create_star_PSF(
            self.star_flux, self.star_coordinates, self.gaussian_std
        )
        header = self.HDR.create_header()

        fits.writeto(
            self.image_dir + self.image_name + ".fits",
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
        self._integrate_spectruns()
        background = self.BGI.create_background_image(self.sky_flux)
        header = self.HDR.create_header()

        fits.writeto(
            self.image_dir + self.image_name + "_BG.fits",
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

        fits.writeto(
            self.image_dir + self.image_name + "_DARK.fits",
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

        fits.writeto(
            self.image_dir + self.image_name + "_BIAS.fits",
            bias,
            overwrite=True,
            header=header,
        )

    def create_radom_image(self, n=10):
        """Create a random star image.

        This function creates an artificial image with a set of random stars.
        The number of star created by the function can be provided to the class. Otherwise,
        the number 10 will be used. Then, all the star images will be sumed with the background image.

        Parameters
        ----------

        n: int, optional
            The number of star in the image.

        Returns
        -------
        Star Image:
            A FITS file with the calculated random artificial image.
        """

        self._integrate_spectruns()
        random_image = self.BGI.create_background_image(self.sky_flux)
        for i in range(n):
            x_coord = randint(50, self.image_size - 50)
            y_coord = randint(1, self.image_size)
            star_flux = randint(10, 4000)
            gaussian_std = randint(1, 8)
            random_image += self.PSF.create_star_PSF(
                star_flux, (x_coord, y_coord), gaussian_std
            )
        header = self.HDR.create_header()

        fits.writeto(
            self.image_dir + self.image_name + ".fits",
            random_image,
            overwrite=True,
            header=header,
        )
