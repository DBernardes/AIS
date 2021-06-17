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
    ccd_operation_mode : dictionary
        Operation mode parameters of the SPARC4 camera

        em_gain : float
            CCD Electron Multiplying gain
        binn : [1, 2]
            Binning of the pixels
        t_exp : float
            Exposure time in seconds
        image_size : int
            Image size in pixels

    ccd_gain : float
        Gain of the CCD in e-/ADU.
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

    def create_star_PSF(
        self,
        star_coordinates,
        gaussian_std,
        ordinary_ray,
        extra_ordinary_ray=0,
    ):
        """Create the star point spread function.

        Parameters
        ----------

        star_coordinates: tuple
            XY star coordinates in the image.

        gaussian_std: int
            Standard deviation of the gaussian 2D distribution in pixels.

        ordinary_ray: float
            Photons/s of the ordinary ray flux.

        extra_ordinary_ray: float, optional
            Photons/s of the extra ordinary ray flux.

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

        gaussian_amplitude = ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
        shape = (image_size, image_size)
        table = Table()
        table["amplitude"] = [gaussian_amplitude]
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std / binn]
        table["y_stddev"] = [gaussian_std / binn]
        table["theta"] = np.radians(np.array([0]))

        self.star_image = make_gaussian_sources_image(shape, table)

        if extra_ordinary_ray != 0:
            gaussian_amplitude = (
                extra_ordinary_ray * t_exp * em_gain * binn ** 2 / ccd_gain
            )
            table["amplitude"] = [gaussian_amplitude]
            table["x_mean"] = [x_coord - 20]
            table["y_mean"] = [y_coord - 20]
            self.star_image += make_gaussian_sources_image(shape, table)

        return self.star_image
