"""
Artificial Images Simulator Package
===================================

The Artificial Images Simulator (AIS) package was developed to generate
artificial star images, similar to those images that would be acquired by
using the acquisition system of the instrument. To accomplish this,
the AIS models as star flux as a 2D gaussian distribution. Then, the star
flux is added to an image with a background level given by counts distribution
of an image of the SPARC4 cameras, as a function of its operation mode.
"""

import read_noise_calc as RNC
from photutils.datasets import make_noise_image
import astropy.io.fits as fits
import numpy as np
import openpyxl

from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image
from sys import exit


__all__ = ['Artificial_Images_Simulator']


class Artificial_Images_Simulator:
    """Create an image cube with the star flux distribution.

    Parameters
    ----------
    star_flux : float
        Number of photons per second of the star
    sky_flux: float
        Number of photons per second of the sky
    gaussian_stddev: int
        Number of pixels of the gaussian standard deviation
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

    ccd_temp: float, optional
        CCD temperature

    serial_number: {9914, 9915, 9916 or 9917}, optional
        CCD serial number

    bias_level: int, optional
        Bias level, in ADU, of the image

    image_dir: str, optional
        Directory where the image should be saved


    Attribute
    ----------
    image_name: str
        Name of the image cube
    dark_current: float
        Dark current in e-/s/pix
    read_noise: float
        Read noise in e-/pix
    gain: float
        CCD pre-amplification gain in e-/ADU
    hdr: list
        Header content of the image


    Yields
    ------
        image cube: int
            An image cube in the FITS format with the star flux distribution

    Notes
    -----
        Explica o código; background; passo-a-passo

    Examples
    --------
        Incluir exemplos

    References
    ----------
    .. [#Bernardes_2018] Bernardes, D. V., Martioli, E., and Rodrigues, C. V., “Characterization of the SPARC4 CCDs”, <i>Publications of the Astronomical Society of the Pacific</i>, vol. 130, no. 991, p. 95002, 2018. doi:10.1088/1538-3873/aacb1e.

    """

    def __init__(self,
                 star_flux,
                 sky_flux,
                 gaussian_stddev,
                 ccd_operation_mode,
                 ccd_temp=-70,
                 serial_number=9914,
                 bias_level=500,
                 image_dir=''):
        """Initialize the class."""
        if type(star_flux) not in [int, float]:
            raise ValueError(f'The star flux must be a number: {star_flux}')
        elif star_flux > 0:
            self.star_flux = star_flux
        else:
            raise ValueError(
                f'The star flux must be greater than zero: {star_flux}')

        if type(sky_flux) not in [int, float]:
            raise ValueError(f'The sky flux must be a number: {sky_flux}')
        elif sky_flux > 0:
            self.sky_flux = sky_flux
        else:
            raise ValueError(
                f'The sky flux must be greater than zero: {sky_flux}')

        if type(gaussian_stddev) is not int:
            raise ValueError(
                f'The gaussian standard deviation must be an integer: {gaussian_stddev}')
        elif gaussian_stddev > 0:
            self.gaussian_stddev = gaussian_stddev
        else:
            raise ValueError(
                f'The gaussian standard deviation must be greater than zero: {gaussian_stddev}')

        self._verify_ccd_operation_mode(ccd_operation_mode)

        if type(ccd_temp) not in [int, float]:
            raise ValueError(
                f'The CCD temperature must be a number: {ccd_temp}')
        elif -70 <= ccd_temp and ccd_temp <= 0:
            self.ccd_temp = ccd_temp
        else:
            raise ValueError(
                f'CCD temperatura out of range [-70 ºC, 0 ºC]: {ccd_temp}')

        if serial_number in [9914, 9915, 9916, 9917]:
            self.serial_number = serial_number
        else:
            raise ValueError(
                f'There is no camera with the provided serial number: {serial_number}')

        if type(bias_level) is not int:
            raise ValueError(
                f'The bias level must be an integer: {bias_level}')
        elif bias_level >= 0:
            self.bias_level = bias_level
        else:
            raise ValueError(f'The bias level must be positive: {bias_level}')

        if type(image_dir) is not str:
            raise ValueError(
                f'The directory path must be a string: {image_dir}')
        else:
            if image_dir != '':
                if '\\' not in image_dir[-1]:
                    image_dir += '\\'
            self.image_dir = image_dir

        self.image_name = ''
        self.dark_current = 0
        self.read_noise = 0
        self.gain = 0
        self.hdr = []

    def _verify_ccd_operation_mode(self, ccd_operation_mode):
        """Verify if the provided CCD operation mode is correct."""
        self.ccd_operation_mode = ccd_operation_mode
        dic_keywords_list = [
            'em_mode', 'em_gain', 'preamp', 'hss', 'bin', 't_exp']
        for key in ccd_operation_mode.keys():
            if key not in dic_keywords_list:
                raise ValueError(
                    f'The name provided is not a CCD parameter: {key}')

        if list(ccd_operation_mode.keys()) != dic_keywords_list:
            raise ValueError(
                'There is a missing parameter of the CCD operation mode')

        if ccd_operation_mode['em_mode'] not in [0, 1]:
            raise ValueError(
                f'Invalid value for the EM mode: {ccd_operation_mode["em_mode"]}')
        if ccd_operation_mode['em_gain'] < 2 or ccd_operation_mode['em_gain'] > 300:
            raise ValueError(
                f'EM gain out of range [2, 300]: {ccd_operation_mode["em_mode"]}')
        if ccd_operation_mode['preamp'] not in [1, 2]:
            raise ValueError(
                f'Invalid value for the pre-amplification: {ccd_operation_mode["preamp"]}')
        if ccd_operation_mode['hss'] not in [0.1, 1, 10, 20, 30]:
            raise ValueError(
                f'Invalid value for the Readout rate: {ccd_operation_mode["hss"]}')
        if ccd_operation_mode['bin'] not in [1, 2]:
            raise ValueError(
                f'Invalid value for the binning: {ccd_operation_mode["bin"]}')
        if ccd_operation_mode['t_exp'] < 1e-5:
            raise ValueError(
                f'Invalid value for the exposure time: {ccd_operation_mode["t_exp"]}')

    def _write_image_mode(self):
        """Write the CCD operation mode to the attributes of the class."""
        self.em_mode = self.ccd_operation_mode['em_mode']
        self.noise_factor = 1
        self.em_gain = 1
        if self.em_mode == 1:
            self.noise_factor = 1.41
            self.em_gain = self.ccd_operation_mode['em_gain']
        self.preamp = self.ccd_operation_mode['preamp']
        self.hss = self.ccd_operation_mode['hss']
        self.bin = self.ccd_operation_mode['bin']
        self.t_exp = self.ccd_operation_mode['t_exp']

    def _configure_gain(self):
        """Configure the CCD gain based on its operation mode."""
        em_mode = self.em_mode
        hss = self.hss
        preamp = self.preamp
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
            r'spreadsheet\Read_noise_and_gain_values.xlsx').active
        self.gain = spreadsheet.cell(tab_index, 5).value

    def _calc_dark_current(self):
        """Calculate the CCD dark current as a function of its temperature.

        The calculation of the DC is based on the model presented in the
        article of the Characterization of the SPARC4 CCDs [#Bernardes_2018]_
        """
        T = self.ccd_temp
        if self.serial_number == 9914:
            self.dark_current = 24.66*np.exp(0.0015*T**2+0.29*T)
        if self.serial_number == 9915:
            self.dark_current = 35.26*np.exp(0.0019*T**2+0.31*T)
        if self.serial_number == 9916:
            self.dark_current = 9.67*np.exp(0.0012*T**2+0.25*T)
        if self.serial_number == 9917:
            self.dark_current = 5.92*np.exp(0.0005*T**2+0.18*T)

    def _calc_read_noise(self):
        """Calculate the read noise the CCD.

        The calculation is performed by providing the CCD operation mode to
        the ReadNoiseCalc package
        """
        RN = RNC.ReadNoiseCalc()
        RN.write_operation_mode(
            self.em_mode, self.em_gain, self.hss, self.preamp, self.bin)
        RN.calc_read_noise()
        self.read_noise = RN.calc_read_noise()

    def _create_image_header(self):
        """Create the image header.

        This functions writes a simple header with the used parameters for
        the CCD operation mode to the image cube file
        """
        hdr = fits.Header()
        hdr['NAXIS1'] = (200, 'length of data axis 1')
        hdr['NAXIS2'] = (200, 'length of data axis 2')
        hdr['EXTEND'] = ('T', 'FITS dataset may contain extensions')
        hdr['COMMENT'] = 'and Astrophysics, volume 376, page 359; bibcode:' + \
            '2001A&A...376..3'
        hdr['ACQMODE'] = ('Single  ', 'Acquisition Mode')
        hdr['READMODE'] = ('Image   ', 'Readout Mode')
        hdr['IMGRECT'] = ('1, 200,200, 1', 'Image Format')
        hdr['HBIN'] = (self.bin, 'Horizontal Binning')
        hdr['VBIN'] = (self.bin, 'Vertical Binning')
        hdr['TRIGGER'] = ('Internal', 'Trigger Mode')
        hdr['EXPOSURE'] = (self.t_exp, 'Total Exposure Time')
        hdr['TEMP'] = (self.ccd_temp, 'Temperature')
        hdr['READTIME'] = (str(1/self.hss)+'E-006', 'Pixel readout time ')
        hdr['VSHIFT'] = ('4.33E-06', 'Vertical Shift Speed')
        hdr['GAIN'] = (self.gain, 'Preamp Gain (e-/ADU)')
        em_mode = 'Conventional'
        if self.em_mode == 1:
            em_mode = 'Electron Multiplying'
        hdr['OUTPTAMP'] = (em_mode, 'Output Amplifier')
        hdr['EMGAIN'] = (self.em_gain, 'Electron Multiplying Gain')
        hdr['PREAMP'] = (str(self.preamp)+'x', 'Pre Amplifier Gain')
        hdr['SERNO'] = (self.serial_number, 'Serial Number')
        hdr['DATE'] = ('2017-07-14T00:00:58',
                       'File Creation Date (YYYY-MM-HHThh:mm:ss)')
        hdr['FRAME'] = ('2017-07-14T00:00:58.642', 'Start of Frame Exposure')
        hdr['IMAGE'] = ('hats-24_I_transito_001', 'Nome do arquivo')
        self.hdr = hdr

    def _configure_image_name(self, include_star_flux=False):
        """Create the image name.

        The image name will be created based on the provided information

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        em_gain = '_G' + str(self.em_gain)
        em_mode = 'CONV'
        if self.em_mode == 1:
            em_mode = 'EM'
        hss = '_HSS' + str(self.hss)
        preamp = '_PA' + str(self.preamp)
        binn = '_B' + str(self.bin)
        t_exp = '_TEXP' + str(self.t_exp)
        self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain

        if include_star_flux:
            star_flux = '_S' + str(self.star_flux)
            self.image_name += star_flux

    def create_image_cube(self):
        """Create the artificial image cube."""
        t_exp = self.t_exp
        em_gain = self.em_gain
        gain = self.gain
        bias = self.bias_level
        dc = self.dark_current*t_exp
        rn = self.read_noise
        sky = self.sky_flux
        nf = self.noise_factor
        binn = self.bin
        gaussian_stddev = self.gaussian_stddev

        gaussian_amplitude = self.star_flux * t_exp * em_gain * binn**2 / gain
        shape = (200, 200)
        table = Table()
        table['amplitude'] = [gaussian_amplitude]
        table['x_mean'] = [100]
        table['y_mean'] = [100]
        table['x_stddev'] = [gaussian_stddev/binn]
        table['y_stddev'] = [gaussian_stddev/binn]
        table['theta'] = np.radians(np.array([0]))

        star_image = make_gaussian_sources_image(shape, table)
        background_level = bias + (dc + sky) * t_exp * em_gain * binn**2 / gain
        image_noise = np.sqrt(rn**2 + (sky + dc)*t_exp *
                              nf**2 * em_gain**2 * binn**2)/gain
        noise_image = make_noise_image(shape,
                                       distribution='gaussian',
                                       mean=background_level,
                                       stddev=image_noise)
        star_image += noise_image
        fits.writeto(self.image_dir + self.image_name + '.fits',
                     star_image, overwrite=True, header=self.hdr)

    def create_bias_image(self):
        """Create a bias image.

        It will be used the same operation mode used for the star image
        """
        gain = self.gain
        bias = self.bias_level
        rn = self.read_noise

        shape = (200, 200)
        background_level = bias
        image_noise = rn
        noise_image = make_noise_image(shape,
                                       distribution='gaussian',
                                       mean=background_level,
                                       stddev=image_noise)/gain
        fits.writeto(self.image_dir + self.image_name + '_BIAS.fits',
                     noise_image, overwrite=True, header=self.hdr)
