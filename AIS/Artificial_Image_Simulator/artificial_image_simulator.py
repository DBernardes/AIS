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
from random import uniform, randint
import astropy.io.fits as fits
import numpy as np


from ..Background_Image import Background_Image
from ..Header import Header
from ..Point_Spread_Function import Point_Spread_Function
from ..Spectral_Response import Atmosphere, Telescope, Channel
from ..Spectral_Energy_Distribution import Source, Sky


__all__ = ['Artificial_Image_Simulator']


class Artificial_Image_Simulator:
    """Artificial Images Simulator class.

    The Artificial Images Simulator (AIS) class was developed to generate artificial star images,
    similar to those images that would be acquired by using the acquisition system of the SPARC4 instrument.
    To accomplish this, the AIS models as star flux as a 2D gaussian distribution.
    Then, the star flux is added to an image with a background level
    given by counts distribution of an image of the SPARC4 cameras, as a function of its operation mode. 

    Example
    --------    

    ccd_operation_mode = {
        'em_mode': 'Conv',
        'em_gain': 1,
        'preamp': 1,
        'readout': 1,
        'binn': 1,
        't_exp': 1,
        'image_size': 100
    }         

    ais = Artificial_Image_Simulator(ccd_operation_mode, channel_id=1, ccd_temperature=-70)

    ais.create_source_sed(calculation_method='blackbody', 
                            magnitude=15, 
                            wavelength_interval=(400, 1100, 1000), 
                            temperature=5700)

    ais.create_sky_sed(moon_phase='new')

    ais.apply_atmosphere_spectral_response()

    ais.apply_telescope_spectral_response()

    ais.apply_sparc4_spectral_response(acquisition_mode='photometry')

    ais.create_artificial_image(image_path='.', star_coordinates=(50,50)) 

    Notes
    -----
        Explicar o código; background; passo-a-passo


    References
    ----------
    .. [#Bernardes_2018] Bernardes, D. V., Martioli, E., and Rodrigues, C. V., “Characterization of the SPARC4 CCDs”, <i>Publications of the Astronomical Society of the Pacific</i>, vol. 130, no. 991, p. 95002, 2018. doi:10.1088/1538-3873/aacb1e.

    """

    def __init__(
            self,
            ccd_operation_mode: dict[str, int | float | str],
            channel_id: int,
            ccd_temperature: float | int) -> None:
        """Initialize the Artificial Image Simulator class.

        Parameters
        ----------
        ccd_operation_mode: dictionary
            A python dictionary with the CCD operation mode. The allowed keywords and values are:            

            * em_mode: [Conv, EM]
                Use the Conv for the Conventional Mode and EM for the Electron Multiplying mode
            * em_gain: float
                Electron Multiplying gain
            * preamp: [1, 2]
                Pre-amplification gain
            * readout: [0.1, 1, 10, 20, 30]
                Readout rate in MHz
            * binn: int
                Number of the binned pixels
            * t_exp: float
                Exposure time in seconds            
            * image_size: int, optional
                Image size in pixels
        channel_id: [1, 2, 3, 4]
            The channel ID. 
        ccd_temperature: float
            The CCD temperature in celsius degrees.        
        """
        self.ccd_operation_mode = ccd_operation_mode
        self.channel_id = channel_id
        self._verify_ccd_operation_mode()
        self.BGI_obj = Background_Image(
            self.ccd_operation_mode, channel_id, ccd_temperature)
        self.SRC_obj = Source()
        self.SKY_obj = Sky()
        self.ATM_obj = Atmosphere()
        self.TEL_obj = Telescope()
        self.CHNNL_obj = Channel(channel_id)
        self.PSF_obj = Point_Spread_Function(
            self.ccd_operation_mode, channel_id)
        self.HDR_obj = Header(self.ccd_operation_mode,
                              ccd_temperature, channel_id)
        return

    @staticmethod
    def _verify_var_in_interval(var, var_name, var_min=0, var_max=2 ** 32):
        if var <= var_min or var >= var_max:
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
                raise ValueError(
                    f"The name provided is not a CCD parameter: {key}")

        keyvalues = list(self.ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError(
                "There is a missing parameter of the CCD operation mode")

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

        self._check_var_in_a_list(
            readout, "readout rate", [0.1, 1, 10, 20, 30])

        self._check_var_in_a_list(preamp, "pre-amplification", [1, 2])

        self._check_var_in_a_list(binn, "binning", [1, 2])

        self._verify_var_in_interval(image_size, "image size", 0, 1025)
        self.image_size = image_size

        self._verify_var_in_interval(t_exp, "exposure time", 1e-5, 84600)
        self.t_exp = t_exp
        return

    @staticmethod
    def _check_var_in_a_list(var, var_name, _list):
        if var not in _list:
            raise ValueError(
                f"The allowed values for the {var_name} are: {_list}")

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
        self.wavelength, self.source_sed = self.SRC_obj.calculate_sed(calculation_method, magnitude,
                                                                      wavelength_interval, temperature, spectral_type)

    def print_available_spectral_types(self) -> None:
        """Print the available spectral types."""
        self.SRC_obj.print_available_spectral_types()
        return

    def create_sky_sed(self, moon_phase: str):
        """Create the Spectral Energy Distribution of the sky.

        Parameters
        ----------
            moon_phase : ['new', 'first quarter', 'third quarter', 'full']
                The phase of the moon.           
        """
        self.sky_sed = self.SKY_obj.calculate_sed(
            moon_phase, self.wavelength)

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

        self.source_sed = self.ATM_obj.apply_spectral_response(
            self.source_sed, self.wavelength, air_mass, sky_condition)

    def apply_telescope_spectral_response(self):
        """Apply the telescope spectral response.

        This functions applies the telescope spectral response on the
        Spectral Energy Distribution of the source and the sky.
        """
        self.source_sed = self.TEL_obj.apply_spectral_response(
            self.source_sed, self.wavelength)
        self.sky_sed = self.TEL_obj.apply_spectral_response(
            self.sky_sed, self.wavelength)

    def apply_sparc4_spectral_response(
        self,
        acquisition_mode: str,
        calibration_wheel: str = "",
        retarder_waveplate: str = "half"
    ):
        """Apply the SPARC4 spectral response.

        This functions applies the SPARC4 spectral response on the
        Spectral Energy Distribution of the source and the sky.

        Parameters
        ----------
        acquisition_mode: ["photometry", "polarimetry"]
            The acquisition mode of the sparc4.

        calibration_wheel: ["polarizer", "depolarizer"], optional
            The position of the calibration wheel.
            This parameter is used only if the acquisition_mode is 'polarimetry'.

        retarder_waveplate: ["half", "quarter"], optional
            The waveplate for polarimetric measurements.
            This parameter is used only if the acquisition_mode is 'polarimetry'.

        """
        self.CHNNL_obj.write_sparc4_operation_mode(
            acquisition_mode, calibration_wheel, retarder_waveplate)
        self.source_sed = self.CHNNL_obj.apply_spectral_response(
            self.source_sed, self.wavelength)
        self.sky_sed = self.CHNNL_obj.apply_spectral_response(
            self.sky_sed, self.wavelength)

    def _integrate_sed(self):
        """Integrate the star and the sky SEDs."""
        self.wavelength *= 1e-9
        self.sky_photons_per_second = np.trapz(self.sky_sed, self.wavelength)
        self.star_photons_per_second = np.trapz(
            self.source_sed, self.wavelength)

    def create_artificial_image(self, image_path: str, star_coordinates: tuple) -> None:
        """
        Create the artificial star image.

        This function creates a FITS file of the artificial star image.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.

        star_coordinates : tuple
            The coordinates in pixels of the star in the image.

        """
        self._integrate_sed()
        ord_ray, extra_ord_ray = self.star_photons_per_second, 0
        if type(self.star_photons_per_second) != np.float64:
            self.sky_photons_per_second = sum(self.sky_photons_per_second)
            ord_ray, extra_ord_ray = self.star_photons_per_second

        background = self.BGI_obj.create_sky_background(
            self.sky_photons_per_second)
        star_psf = self.PSF_obj.create_star_image(
            star_coordinates,
            ord_ray, extra_ord_ray)
        self._create_image_name(image_path)
        header = self.HDR_obj.create_header()
        header['OBSTYPE'] = 'OBJECT'
        header['FILENAME'] = self.image_name
        header['SHUTTER'] = 'OPEN'

        file = os.path.join(image_path, self.image_name)

        fits.writeto(
            file,
            background + star_psf,
            header=header,
        )

    def create_background_image(self, image_path: str, images: int = 1):
        """Create the background image.

        This function creates the background image, similar to those
        acquired by the SPARC4 cameras.  

        Parameters  
        ----------

        image_path : str
            The path where the FITS file will be saved.

        images : int, optional
            The number of images to be created. The default is 1.
        """

        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            self._integrate_sed()
            background = self.BGI_obj.create_sky_background(
                self.sky_photons_per_second)
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header['FILENAME'] = self.image_name
            header['OBSTYPE'] = 'FLAT'  # ?
            header['SHUTTER'] = 'OPEN'

            file = os.path.join(image_path, self.image_name)

            fits.writeto(
                file,
                background,
                header=header
            )

    def create_dark_image(self, image_path: str, images: int = 1):
        """
        Create a dark image.

        This function creates a dark image, similar to those
        acquired by the SPARC4 cameras.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.

        images : int, optional
            The number of dark images to be created. The default is 1.

        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            dark_image = self.BGI_obj.create_dark_background()
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header['FILENAME'] = self.image_name
            header['OBSTYPE'] = 'DARK'
            image_name = os.path.join(
                image_path, self.image_name)

            fits.writeto(
                image_name,
                dark_image,
                header=header,
            )

    def create_bias_image(self, image_path: str, images: int = 1):
        """
        Create a bias image.

        This function creates a bias image, similar to those
        acquired by the SPARC4 cameras.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.       

        images : int, optional
            The number of bias images to be created. The default is 1.
        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            bias = self.BGI_obj.create_bias_background()
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header["EXPTIME"] = 1e-5
            header['FILENAME'] = self.image_name
            header['OBSTYPE'] = 'ZERO'

            image_name = os.path.join(image_path, self.image_name)

            fits.writeto(
                image_name,
                bias,
                header=header,
            )

    def create_random_image(self, image_path, number_stars=10):
        """
        Create a random star image.

        This function creates an artificial image with a set of random stars.
        The number of star created by the function can be provided to the class. Otherwise,
        the number 10 will be used. Then, all the star images will be sumed with the background image.

        Parameters
        ----------

        image_path : str
            The path where the FITS file will be saved.

        number_stars: int, optional
            The number of stars in the image.

        Returns
        -------
        Star Image:
            A FITS file with the calculated random artificial image.
        """
        self._integrate_sed()
        ord_ray, extra_ord_ray = self.star_photons_per_second, 0
        if type(self.star_photons_per_second) != np.float64:
            self.sky_photons_per_second = sum(self.sky_photons_per_second)
            ord_ray, extra_ord_ray = self.star_photons_per_second

        image_size = self.ccd_operation_mode['image_size']
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
            )

        background = self.BGI_obj.create_sky_background(
            self.sky_photons_per_second)
        self._create_image_name(image_path)
        header = self.HDR_obj.create_header()
        header['OBSTYPE'] = 'OBJECT'
        header['FILENAME'] = self.image_name
        header['SHUTTER'] = 'OPEN'
        file = os.path.join(image_path, self.image_name)
        fits.writeto(
            file,
            background + star_image,
            header=header,
        )

    def create_flat_image(self, image_path: str, images: int = 1) -> None:
        """
        Create a flat image.

        This function creates a flat image, similar to those
        acquired by the SPARC4 cameras.

        Parameters
        ----------
        image_path : str
            The path where the FITS file will be saved.    

        images : int, optional
            The number of flat images to be created. The default is 1.   
        """
        if images < 1:
            raise ValueError("The number of images must be greater than 0.")

        for _ in range(images):
            flat_image = self.BGI_obj.create_flat_background()
            header = self.HDR_obj.create_header()
            self._create_image_name(image_path)
            header['FILENAME'] = self.image_name
            header['OBSTYPE'] = 'FLAT'
            header['SHUTTER'] = 'OPEN'

            image_name = os.path.join(
                image_path, self.image_name)

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
            f"%Y%m%d_s4c{self.channel_id}_{index + 1:06}.fits")
        return
