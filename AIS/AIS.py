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

from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image


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

    serial_number: int, optional
        CCD serial number

    bias_level: int, optional
        Bias level, in ADU, of the image

    image_dir: str, optional
        Directory where the image should be saved


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

    def __init__(self, star_flux, sky_flux,
                 gaussian_stddev,
                 ccd_operation_mode,
                 ccd_temp=-70,
                 serial_number=9916,
                 bias_level=500,
                 image_dir=''):
        """Initialize the class."""
        # star flux in photons/s
        self.star_flux = star_flux
        # sky flux in photons/s
        self.sky_flux = sky_flux
        # standard deviaion of the 2D-Gaussian Distribution in pixels
        self.gaussian_stddev = gaussian_stddev
        # operation mode of the CCD.
        self.ccd_operation_mode = ccd_operation_mode
        # CCD temperature. It should be between 0 ºC to -70 ºC
        self.ccd_temp = ccd_temp
        # Serial number of the CCD. For the SPARC4 cameras, they would be
        # 9914, 9915, 9916, or 9917.
        self.serial_number = serial_number
        # Bias level in analogical-to digital units for the pixels of the
        # created image
        self.bias_level = bias_level
        # Directory where the image should be saved
        if image_dir != '':
            if '\\' not in image_dir[-1]:
                image_dir += '\\'
        self.image_dir = image_dir
        # Name of the createde image. It is automatically generated.
        self.image_name = ''

        # Dark current of the CCD.
        self.dark_current = 0
        # Read noise of the CCD.
        self.read_noise = 0
        # Gain of the CCD in e-/ADU
        self.gain = 0
        # Image header
        self.hdr = []

    def write_image_mode(self):
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

        # Calculated the gain, dark current and the read noise of the CCD for
        # the provided operation mode
        self.set_gain()
        self.set_dc()
        self.calc_RN()

    def set_gain(self):
        """Set the CCD gain based on its operation mode."""
        em_mode = self.em_mode
        hss = self.hss
        preamp = self.preamp
        gain = 0
        if em_mode == 1:
            if hss == 30:
                if preamp == 1:
                    gain = 17.2
                if preamp == 2:
                    gain = 5.27
            if hss == 20:
                if preamp == 1:
                    gain = 16.4
                if preamp == 2:
                    gain = 4.39
            if hss == 10:
                if preamp == 1:
                    gain = 16.0
                if preamp == 2:
                    gain = 3.96
            if hss == 1:
                if preamp == 1:
                    gain = 15.9
                if preamp == 2:
                    gain = 3.88
        else:
            if hss == 1:
                if preamp == 1:
                    gain = 3.37
                if preamp == 2:
                    gain = 0.8
            if hss == 0.1:
                if preamp == 1:
                    gain = 3.35
                if preamp == 2:
                    gain = 0.8
        self.gain = gain

    def calc_dark_current(self):
        """Calculate the CCD dark current.

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

    def calc_RN(self):
        """Calculate the read noise the CCD.

        The calculation is performed by providing the CCD operation mode to
        the ReadNoiseCalc package
        """
        RN = RNC.ReadNoiseCalc()
        RN.write_operation_mode(
            self.em_mode, self.em_gain, self.hss, self.preamp, self.bin)
        RN.calc_read_noise()
        self.read_noise = RN.calc_read_noise()

    def create_image_header(self):
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

    def write_image_name(self, include_star_flux=False):
        """Create the image name.

        It will be created the image name based on the provided CCD operation
        mode

        Parameters
        ----------
        include_star_flux: bool, optional
            Indicate if it is needed to include the star flux value in the
            image name
        """
        em_gain = '_G' + str(self.em_gain)
        em_mode = 'CONV_'
        if self.em_mode == 1:
            em_mode = 'EM_'
        hss = str(self.hss) + 'MHz'
        preamp = '_PA' + str(self.preamp)
        binn = '_B' + str(self.bin)
        t_exp = '_TEXP' + str(self.t_exp)
        self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain

        if include_star_flux:
            star_flux = '_S' + str(self.star_flux)
            self.image_name += star_flux

    def create_artificial_image(self):
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

        Create a bias image with the same operation mode used for the star
        image
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
