"""
Background Image 
================

This is the Background Image Class used to generate a back ground image like
a bias image acquired by the SPARC4 cameras.
"""

from FC import Flux_Calculation
from Telescope_SR import Telescope_Spectral_Response
from Atmosphere_SR import Atmosphere_Spectral_Response
from CHC import Abstract_Channel_Creator
import numpy as np


from photutils.datasets import make_noise_image
from astropy.table import Table


class Background_Image:
    """Background Image Class.

    Parameters
    ----------
    Abstract_Channel_Creator : object
        An object of the Channel Creator Class. This object is used to
        calculate the contribution of the SPARC4 instrument in the sky flux

    ccd_operation_mode: dictionary
        A python dictionary with the CCD operation mode.
        The allowed keywords values for the dictionary are

        * em_mode: {0, 1}

           Use the 0 for the Conventional Mode and 1 for the EM Mode

        * preamp: {1, 2}

           Pre-amplification

        * hss: {0.1, 1, 10, 20, 30}

           Horizontal Shift Speed (readout rate) in MHz

        * bin: int

           Number of the binned pixels

    ccd_gain : float
        CCD gain in e-/ADU.
        """

    def __init__(self, Abstract_Channel_Creator,
                 ccd_operation_mode,
                 ccd_gain):
        """Initialize the Background Image class."""
        self.CHC = Abstract_Channel_Creator
        self.FC = Flux_Calculation()
        self.TSR = Telescope_Spectral_Response()
        self.ASR = Atmosphere_Spectral_Response()
        self.ccd_gain = ccd_gain

        self.em_gain = ccd_operation_mode['em_gain']
        self.binn = ccd_operation_mode['binn']
        self.t_exp = ccd_operation_mode['t_exp']
        self.preamp = ccd_operation_mode['preamp']
        self.hss = ccd_operation_mode['hss']

        self.BIAS_LEVEL = 100
        self.NOISE_FACTOR = 1.0
        if ccd_operation_mode['em_mode'] == 1:
            self.NOISE_FACTOR = 1.4

        self._calculate_sky_flux()
        self._calculate_dark_current()
        self._calculate_read_noise(ccd_operation_mode)

    def _calculate_sky_flux(self):
        self.sky_flux = self.FC.calc_sky_flux()

    def _calculate_dark_current(self):
        self.dark_current = self.CHC.calc_dark_current()

    def _calculate_read_noise(self, ccd_operation_mode):
        self.read_noise = self.CHC.calculate_read_noise(ccd_operation_mode)

    def create_background_image(self):
        """Create the background image.

        This functions creates a background image with a background level given
        by the ccd operation mode, the sky flux, the dc noise, and the bias
        level. Over this image there is a noise given by a gaussian
        distribution over the read noise, dc noise, and sky noise. Also, the
        extra noise of the EM amplification is considered.

        Returns
        -------
        noise_image: array like
            A background image for the respective CCD operation mode, the sky
            flux, and the dc level.

        """
        t_exp = self.t_exp
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        bias = self.BIAS_LEVEL
        dc = self.dark_current * t_exp
        rn = self.read_noise
        sky = self.sky_flux
        nf = self.NOISE_FACTOR
        binn = self.binn

        shape = (200, 200)
        table = Table()
        table['amplitude'] = [0]
        table['x_mean'] = [100]
        table['y_mean'] = [100]
        table['x_stddev'] = [0]
        table['y_stddev'] = [0]
        table['theta'] = [0]

        background_level = bias \
            + (dc + sky) * t_exp * em_gain * binn**2 / ccd_gain

        noise = np.sqrt(rn**2 + (sky + dc)*t_exp
                        * nf**2 * em_gain**2 * binn**2)/ccd_gain

        self.background_image = make_noise_image(shape,
                                                 distribution='gaussian',
                                                 mean=background_level,
                                                 stddev=noise)

        return self.background_image
