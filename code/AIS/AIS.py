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

from sys import exit
import openpyxl
import astropy.io.fits as fits
from PSF import Point_Spread_Function
from BGI import Background_Image
from HDR import Header
from CHC import (Concrete_Channel_1,
                 Concrete_Channel_2,
                 Concrete_Channel_3,
                 Concrete_Channel_4)


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

    channel: {1, 2, 3, 4}, optional
        The SPARC4 channel

    gaussian_stddev: int
        Number of pixels of the gaussian standard deviation

    star_position: list
        XY coordinates of the position of the star

    image_size: int, optional
        Image size in pixels

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

    def __init__(self,
                 star_magnitude,
                 sky_magnitude,
                 ccd_operation_mode,
                 channel=1,
                 gaussian_std=3,
                 star_position=[100, 100],
                 image_size=200,
                 bias_level=500,
                 sparc4_operation_mode='phot',
                 image_dir=''):
        """Initialize the class."""
        if type(star_magnitude) not in [int, float]:
            raise ValueError('The star flux must be a number: '
                             + f'{star_magnitude}')
        elif star_magnitude <= 0:
            raise ValueError(
                f'The star flux must be greater than zero: {star_magnitude}')
        else:
            self.star_magnitude = star_magnitude

        if type(sky_magnitude) not in [int, float]:
            raise ValueError(f'The sky flux must be a number: {sky_magnitude}')
        elif sky_magnitude <= 0:
            raise ValueError(
                f'The sky flux must be greater than zero: {sky_magnitude}')
        else:
            self.sky_magnitude = sky_magnitude

        if channel in [1, 2, 3, 4]:
            self.channel = channel
        else:
            raise ValueError(
                'There is no camera with the provided'
                + f'serial number: {channel}')

        if type(gaussian_std) is not int:
            raise ValueError(
                'The gaussian standard deviation must be'
                + f'an integer: {gaussian_std}')
        elif gaussian_std <= 0:
            raise ValueError(
                r'The gaussian standard deviation must be greater'
                + f'than zero: {gaussian_std}')
        else:
            self.gaussian_std = gaussian_std

        for coord in star_position:
            if type(coord) is not int:
                raise ValueError(
                    'The star coordinates must be an integer: {coord}')
            elif coord <= 0:
                raise ValueError(
                    'The star coordinates must be greater than zero: {coord}')
            else:
                self.star_position = star_position

        if type(image_size) is not int:
            raise ValueError(
                f'The image size must be an integer: {image_size}')
        elif image_size <= 0:
            raise ValueError(
                f'The image size must be greater than zero: {image_size}')
        else:
            self.image_size = image_size

        if type(bias_level) is not int:
            raise ValueError(
                f'The bias level must be an integer: {bias_level}')
        elif bias_level <= 0:
            raise ValueError(f'The bias level must be positive: {bias_level}')
        else:
            self.bias_level = bias_level

        if type(sparc4_operation_mode) is not str:
            raise ValueError(
                r'The SPARC4 operation mode must be a string: '
                + f'{sparc4_operation_mode}')
        elif sparc4_operation_mode not in ['phot', 'pol']:
            raise ValueError('The SPARC4 operation mode must be "phot" '
                             + f'or "pol": {sparc4_operation_mode}')
        else:
            self.sparc4_operation_mode = sparc4_operation_mode

        if type(image_dir) is not str:
            raise ValueError(
                f'The directory path must be a string: {image_dir}')
        else:
            if image_dir != '':
                if '/' not in image_dir[-1]:
                    image_dir += '/'
            self.image_dir = image_dir

        self._verify_ccd_operation_mode(ccd_operation_mode)
        self._configure_gain(ccd_operation_mode)
        self._configure_image_name(ccd_operation_mode)

        CHC = 0
        ccd_temp = ccd_operation_mode['ccd_temp']
        if channel == 1:
            CHC = Concrete_Channel_1(ccd_temp,
                                     sparc4_acquisition_mode='phot')
        elif channel == 2:
            CHC = Concrete_Channel_2(ccd_temp,
                                     sparc4_acquisition_mode='phot')
        elif channel == 3:
            CHC = Concrete_Channel_3(ccd_temp,
                                     sparc4_acquisition_mode='phot')
        elif channel == 4:
            CHC = Concrete_Channel_4(ccd_temp,
                                     sparc4_acquisition_mode='phot')
        self.CHC = CHC
        self.PSF = Point_Spread_Function(
            CHC, ccd_operation_mode, self.ccd_gain, self.gaussian_std)
        self.BGI = Background_Image(CHC, ccd_operation_mode, self.ccd_gain,
                                    self.bias_level)
        self.HDR = Header(ccd_operation_mode, self.ccd_gain,
                          CHC.get_serial_number())

    def _verify_ccd_operation_mode(self, ccd_operation_mode):
        """Verify if the provided CCD operation mode is correct."""
        dic_keywords_list = [
            'binn', 'ccd_temp', 'em_gain', 'em_mode', 'hss', 'preamp', 't_exp']

        for key in ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(
                    f'The name provided is not a CCD parameter: {key}')

        keyvalues = list(ccd_operation_mode.keys())
        keyvalues.sort()
        if keyvalues != dic_keywords_list:
            raise ValueError(
                'There is a missing parameter of the CCD operation mode')

        em_mode = ccd_operation_mode['em_mode']
        em_gain = ccd_operation_mode['em_gain']
        hss = ccd_operation_mode['hss']
        preamp = ccd_operation_mode['preamp']
        binn = ccd_operation_mode['binn']
        t_exp = ccd_operation_mode['t_exp']
        ccd_temp = ccd_operation_mode['ccd_temp']

        if em_mode not in [0, 1]:
            raise ValueError(
                f'Invalid value for the EM mode: {em_mode}')
        if em_mode == 0:
            if em_gain != 1:
                raise ValueError(
                    'The EM Gain must be 1 for the Conventional'
                    + f' Mode: {em_gain}')
        else:
            if em_gain not in [float, int]:
                raise ValueError(
                    f'The EM gain must be a number: {em_gain}')
            elif em_gain < 2 or em_gain > 300:
                raise ValueError(
                    f'EM gain out of range [2, 300]: {em_gain}')

        if preamp not in [1, 2]:
            raise ValueError(
                f'Invalid value for the pre-amplification: {preamp}')

        if hss not in [0.1, 1, 10, 20, 30]:
            raise ValueError(
                f'Invalid value for the Readout rate: {hss}')

        if binn not in [1, 2]:
            raise ValueError(
                f'Invalid value for the binning: {bin}')

        if type(t_exp) not in [float, int]:
            raise ValueError(
                f'The exposure time must be a number: {t_exp}')
        elif ccd_operation_mode['t_exp'] < 1e-5:
            raise ValueError(
                f'Invalid value for the exposure time: {t_exp}')

        if type(ccd_temp) not in [float, int]:
            raise ValueError(
                f'The CCD temperature must be a number: {ccd_temp}')
        if ccd_temp < -80 or ccd_temp > 20:
            raise ValueError(
                f'CCD temperature out of range [-80, 20]: {ccd_temp}')

    def get_channel_ID(self):
        """Return the ID for the respective SPARC4 channel."""
        return self.CHC.get_channel_ID()

    def _configure_image_name(self, ccd_operation_mode,
                              include_star_mag=False):
        """Create the image name.

        The image name will be created based on the provided information

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        dic = ccd_operation_mode
        em_gain = '_G' + str(dic['em_gain'])
        em_mode = 'CONV'
        if dic['em_mode'] == 1:
            em_mode = 'EM'
        hss = '_HSS' + str(dic['hss'])
        preamp = '_PA' + str(dic['preamp'])
        binn = '_B' + str(dic['binn'])
        t_exp = '_TEXP' + str(dic['t_exp'])
        self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain

        if include_star_mag:
            star_flux = '_S' + str(self.star_magnitude)
            self.image_name += star_flux

    def _configure_gain(self, ccd_operation_mode):
        """Configure the CCD gain based on its operation mode."""
        em_mode = ccd_operation_mode['em_mode']
        hss = ccd_operation_mode['hss']
        preamp = ccd_operation_mode['preamp']
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
            raise ValueError('Unexpected value for the readout rate: {hss}')
        if preamp == 2:
            tab_index += 2

        spreadsheet = openpyxl.load_workbook(
            f'./RNC/spreadsheet/Channel {self.channel}'
            + '/Read_noise_and_gain_values.xlsx').active
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
        background = self.BGI.create_background_image()
        star_PSF = self.PSF.create_star_PSF()
        header = self.HDR.create_header()

        fits.writeto(self.image_dir + self.image_name + '.fits',
                     background + star_PSF, overwrite=True, header=header)
