"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class calculates the star flux distribution
based on a gaussian 2D distributiopn
"""

from FC import Flux_Calculation
from Telescope_SR import Telescope_Spectral_Response
from Atmosphere_SR import Atmosphere_Spectral_Response
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image
import numpy as np


class Point_Spread_Function:
    """Point Spread Function Class.

    This class creates the image of the star. The PSF of the
    object is based on a 2D gaussian distribution. Initially, the star flux is
    calculated using the library Flux_Calculation. Over this flux, the
    contribution of the atmosphere, the telescope, and the SPARC4 are
    considered as a function of their spectral response.
    """

    def __init__(self,
                 Abstract_Channel_Creator,
                 ccd_operation_mode,
                 ccd_gain,
                 gaussian_std):
        """Initialize the class.

        Parameters
        ----------
        Abstract_Channel_Creator : object
            A chield object of the Channel Creator class. This object is used
            to calculated the contribution of the instrument spectral response
            over the star flux as a function of the instrument channel.
        ccd_gain : float
            Gain of the CCD in e-/ADU.

        Returns
        -------
        None.

        """
        self.CHC = Abstract_Channel_Creator
        self.FC = Flux_Calculation()
        self.TSR = Telescope_Spectral_Response()
        self.ASR = Atmosphere_Spectral_Response()

        self.em_gain = ccd_operation_mode['em_gain']
        self.binn = ccd_operation_mode['binn']
        self.t_exp = ccd_operation_mode['t_exp']
        self.ccd_gain = ccd_gain
        self.gaussian_std = gaussian_std

        self._calculate_star_flux()

    def _calculate_star_flux(self):
        self.star_flux = self.FC.calc_star_flux()

    def calculate_star_PSF(self):
        """Create the artificial image cube."""
        t_exp = self.t_exp
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        binn = self.binn
        gaussian_std = self.gaussian_std

        gaussian_amplitude = self.star_flux \
            * t_exp * em_gain * binn**2 / ccd_gain
        shape = (200, 200)
        table = Table()
        table['amplitude'] = [gaussian_amplitude]
        table['x_mean'] = [100]
        table['y_mean'] = [100]
        table['x_stddev'] = [gaussian_std/binn]
        table['y_stddev'] = [gaussian_std/binn]
        table['theta'] = np.radians(np.array([0]))

        self.star_image = make_gaussian_sources_image(shape, table)

        return self.star_image
