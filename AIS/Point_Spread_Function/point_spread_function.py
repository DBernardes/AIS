"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class calculates the star flux distribution
based on a gaussian 2D distributiopn
"""


import numpy as np
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image


class Point_Spread_Function:
    """Point Spread Function Class.

    This class creates the image of the star. The PSF of the
    object is based on a 2D gaussian distribution. Initially, the star flux is
    calculated using the library Flux_Calculation. Over this flux, the
    contribution of the atmosphere, the telescope, and the SPARC4 are
    considered as a function of their spectral response.


    Parameters
    ----------
    Abstract_Channel_Creator : object
        A chield object of the Channel Creator class. This object is used
        to calculated the contribution of the instrument spectral response
        over the star flux as a function of the instrument channel.
    ccd_gain : float
        Gain of the CCD in e-/ADU.
    gaussian_std: int
        Standard deviation of the gaussian 2D distribution.

    Returns
    -------
    None.
    """

    def __init__(
        self,
        Abstract_Channel_Creator,
        ccd_operation_mode,
        ccd_gain,
    ):
        """Initialize the class."""
        self.CHC = Abstract_Channel_Creator

        self.em_gain = ccd_operation_mode["em_gain"]
        self.binn = ccd_operation_mode["binn"]
        self.t_exp = ccd_operation_mode["t_exp"]
        self.image_size = ccd_operation_mode["image_size"]
        self.ccd_gain = ccd_gain

    def create_star_PSF(self, star_flux, star_coordinates, gaussian_std):
        """Create the artificial image cube.

        Parameters
        ----------

        star_flux: float
            Photons/s radiated by the star.

        star_coordinates: tuple
            XY star coordinates in the image.

        gaussian_std: int
            Standard deviation of the gaussian 2D distribution.

        Returns
        -------

        star_image: array like
            The Point Spred Function of the star for the respective CCD operation mode.
        """
        t_exp = self.t_exp
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        binn = self.binn
        image_size = self.image_size
        x_coord = star_coordinates[0]
        y_coord = star_coordinates[1]

        gaussian_amplitude = star_flux * t_exp * em_gain * binn ** 2 / ccd_gain
        shape = (image_size, image_size)
        table = Table()
        table["amplitude"] = [gaussian_amplitude]
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std / binn]
        table["y_stddev"] = [gaussian_std / binn]
        table["theta"] = np.radians(np.array([0]))

        self.star_image = make_gaussian_sources_image(shape, table)

        return self.star_image
