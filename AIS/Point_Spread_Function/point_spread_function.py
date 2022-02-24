"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class calculates the star flux distribution
based on a gaussian 2D distributiopn
"""


import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image, make_noise_image


class Point_Spread_Function:
    """
    Point Spread Function Class.

    This class creates the image of the star. The PSF of the
    object is based on a 2D gaussian distribution. Initially, the star flux is
    calculated using the library Flux_Calculation. Over this flux, the
    contribution of the atmosphere, the telescope, and the SPARC4 are
    considered as a function of their spectral response.


    Parameters
    ----------
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

    _SPARC4_POL_SEPARATION = 40  # pix
    _SPARC4_PLATE_SCALE = 0.35  # arcsec/pix

    def __init__(self, ccd_operation_mode, ccd_gain, seeing):
        """Initialize the class."""
        self.em_gain = ccd_operation_mode["em_gain"]
        self.binn = ccd_operation_mode["binn"]
        self.t_exp = ccd_operation_mode["t_exp"]
        self.image_size = ccd_operation_mode["image_size"]
        self.ccd_gain = ccd_gain
        self.seeing = seeing

        self.star_image = []

    @staticmethod
    def _make_noise_image(shape, table, amplitude):
        table["amplitude"] = [1]
        noise_image = (
            make_noise_image(shape, distribution="poisson", mean=amplitude) - amplitude
        )
        noise_image *= make_gaussian_sources_image(shape, table)

        return noise_image

    def create_star_psf(
        self,
        star_coordinates,
        ordinary_ray,
        extra_ordinary_ray=0,
    ):
        """
        Create the star point spread function.

        Parameters
        ----------

        star_coordinates: tuple
            XY star coordinates in the image.

        ordinary_ray: float
            Photons/s of the ordinary ray flux.

        extra_ordinary_ray: float, optional
            Photons/s of the extra ordinary ray flux.

        Returns
        -------

        star_image: array like
            The Point Spred Function of the star for the respective CCD operation mode.
        """

        x_coord = star_coordinates[0]
        y_coord = star_coordinates[1]
        gaussian_std = self.seeing / self._SPARC4_PLATE_SCALE

        gaussian_amplitude = (
            ordinary_ray * self.t_exp * self.em_gain * self.binn ** 2 / self.ccd_gain
        )
        shape = (self.image_size, self.image_size)

        table = Table()
        table["amplitude"] = [gaussian_amplitude]
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std / self.binn]
        table["y_stddev"] = [gaussian_std / self.binn]
        table["theta"] = np.radians(np.array([0]))

        star_image = make_gaussian_sources_image(shape, table)
        star_image += self._make_noise_image(shape, table, gaussian_amplitude)

        if extra_ordinary_ray != 0:
            gaussian_amplitude = (
                extra_ordinary_ray
                * self.t_exp
                * self.em_gain
                * self.binn ** 2
                / self.ccd_gain
            )
            table["amplitude"] = [gaussian_amplitude]
            table["x_mean"] = [x_coord - self._SPARC4_POL_SEPARATION]
            table["y_mean"] = [y_coord - self._SPARC4_POL_SEPARATION]
            star_image += make_gaussian_sources_image(shape, table)
            star_image += self._make_noise_image(shape, table, gaussian_amplitude)

        return star_image
