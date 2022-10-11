"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class models the star flux distribution
as a gaussian 2D distribution
"""

from array import array
from functools import cached_property
import os
import pandas as pd
import numpy as np
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image, make_noise_image


class Point_Spread_Function:
    """
    Point Spread Function Class.

    This class creates the image of the star flux distribution. The PSF of the
    object is based on a 2D gaussian distribution. Initially, the star flux is
    calculated using the library Flux_Calculation. Over this flux, the
    contribution of the atmosphere, the telescope, and the SPARC4 are
    considered as a function of their spectral response.


    Parameters
    ----------
    ccd_operation_mode : dictionary
        Operation mode parameters of the SPARC4 camera

        em_mode: ['EM', 'Conv']
            The Electron Multiplying mode of the camera.
        em_gain: float
            CCD Electron Multiplying gain.
        readout: [30, 20, 10, 1, 0.1]
            The readout rate of the pixels in MHz.
        preamp:
            The pre-amplification gain.
        binn: [1, 2]
            Binning of the pixels
        t_exp: float
            Exposure time in seconds
        image_size : int
            Image size in pixels   
    seeing: float
        The size of the seeing disk in arcsec.
    channel: int
        The channel related to the camera.
    """

    _SPARC4_POL_SEPARATION = 20  # pix
    _SPARC4_PLATE_SCALE = 0.35  # arcsec/pix
    _SPREADSHEET_PATH = os.path.join(
        'AIS', 'Point_Spread_Function', 'preamp_gains.csv')

    def __init__(self, ccd_operation_mode: dict, channel: int):
        """Initialize the class."""
        self.ccd_operation_mode = ccd_operation_mode
        self.channel = channel
        self.get_ccd_gain()

    def get_ccd_gain(self):
        idx_tab = 0
        readout = self.ccd_operation_mode['readout']
        if self.ccd_operation_mode['em_mode'] == 'EM':
            idx_tab = [30, 20, 10, 1].index(readout) * 2
        else:
            idx_tab = [1, 0.1].index(readout) * 2 + 8
        idx_tab += self.ccd_operation_mode['preamp'] - 1
        ss = pd.read_csv(self._SPREADSHEET_PATH)
        self.ccd_gain = ss[f'CH{self.channel}'][idx_tab]

    @staticmethod
    def _make_noise_image(table, shape):
        amplitude = table["amplitude"]
        table["amplitude"] = [1]
        noise_image = (
            make_noise_image(shape, "poisson", amplitude) - amplitude
        )
        noise_image *= make_gaussian_sources_image(shape, table)

        return noise_image

    def _create_table(self, star_coordinates, seeing):
        table = Table()
        binn = self.ccd_operation_mode['binn']
        gaussian_std = seeing / self._SPARC4_PLATE_SCALE
        x_coord = star_coordinates[0]
        y_coord = star_coordinates[1]
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std / binn]
        table["y_stddev"] = [gaussian_std / binn]
        table["theta"] = np.radians(np.array([0]))
        self.table = table

    def create_star_image(
        self,
        star_coordinates: tuple,
        ordinary_ray: float,
        extra_ordinary_ray: float = 0,
        seeing: float = 1
    ):
        """Create the star image.

        Parameters
        ----------

        star_coordinates: tuple
            XY star coordinates in the image.

        ordinary_ray: float
            Photons/s of the ordinary ray.

        extra_ordinary_ray: float, optional
            Photons/s of the extra ordinary ray.

        seeing: float, optional
            The size of the seeing disk.

        Returns
        -------

        star_image: array like
            The Point Spred Function of the star for the respective CCD operation mode.
        """
        self._create_table(star_coordinates, seeing)
        star_image = self._create_image_ordinary_ray(
            ordinary_ray)
        if extra_ordinary_ray != 0:
            star_image += self._create_image_extra_ordinary_ray(
                extra_ordinary_ray)
        return star_image

    def _create_image_ordinary_ray(self, ordinary_ray):
        t_exp = self.ccd_operation_mode['t_exp']
        em_gain = self.ccd_operation_mode['em_gain']
        binn = self.ccd_operation_mode['binn']
        image_size = self.ccd_operation_mode['image_size']
        gaussian_amplitude = (
            ordinary_ray * t_exp * em_gain * binn ** 2 / self.ccd_gain
        )
        self.table["amplitude"] = [gaussian_amplitude]

        shape = (image_size, image_size)
        star_image = make_gaussian_sources_image(shape, self.table)
        star_image += self._make_noise_image(self.table, shape)
        return star_image

    def _create_image_extra_ordinary_ray(self, extra_ordinary_ray):
        t_exp = self.ccd_operation_mode['t_exp']
        em_gain = self.ccd_operation_mode['em_gain']
        binn = self.ccd_operation_mode['binn']
        image_size = self.ccd_operation_mode['image_size']
        gaussian_amplitude = (
            extra_ordinary_ray
            * t_exp
            * em_gain
            * binn ** 2
            / self.ccd_gain
        )
        self.table["amplitude"] = [gaussian_amplitude]
        self.table["x_mean"] -= self._SPARC4_POL_SEPARATION
        self.table["y_mean"] -= self._SPARC4_POL_SEPARATION
        shape = (image_size, image_size)
        star_image = make_gaussian_sources_image(shape, self.table)
        star_image += self._make_noise_image(self.table, shape)

        return star_image
