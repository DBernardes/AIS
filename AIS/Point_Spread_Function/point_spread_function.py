"""
Point Spread Function Class
===========================

The Point Spread Function (PSF) class models the star flux distribution
as a gaussian 2D distribution
"""

import os
from array import array
from functools import cached_property
from math import pi

import numpy as np
import pandas as pd
from astropy.modeling.models import Gaussian2D
from astropy.table import Table
from numpy import ndarray
from photutils.datasets import make_model_image, make_noise_image

__all__ = ["Point_Spread_Function"]


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

    SPARC4_POL_SEPARATION = 20  # pix
    SPARC4_PLATE_SCALE = 0.35  # arcsec/pix
    SPREADSHEET_PATH = os.path.join("AIS", "Point_Spread_Function", "preamp_gains.csv")

    def __init__(self, ccd_operation_mode: dict, channel: int) -> None:
        """Initialize the class."""
        self.ccd_operation_mode = ccd_operation_mode
        image_size = self.ccd_operation_mode["image_size"]
        self.shape = (image_size, image_size)
        self.channel = channel
        self.get_ccd_gain()

    def get_ccd_gain(self) -> None:
        idx_tab = 0
        readout = self.ccd_operation_mode["readout"]
        if self.ccd_operation_mode["em_mode"] == "EM":
            idx_tab = [30, 20, 10, 1].index(readout) * 2
        else:
            idx_tab = [1, 0.1].index(readout) * 2 + 8
        idx_tab += self.ccd_operation_mode["preamp"] - 1
        ss = pd.read_csv(self.SPREADSHEET_PATH)
        self.ccd_gain = ss[f"CH{self.channel}"][idx_tab]

    def _create_table(self, star_coordinates, seeing) -> Table:
        table = Table()
        binn = self.ccd_operation_mode["binn"]
        gaussian_std = seeing / (self.SPARC4_PLATE_SCALE * binn * 2)
        x_coord = star_coordinates[0]
        y_coord = star_coordinates[1]
        table["x_mean"] = [x_coord]
        table["y_mean"] = [y_coord]
        table["x_stddev"] = [gaussian_std]
        table["y_stddev"] = [gaussian_std]
        table["theta"] = np.radians(np.array([0]))
        self.table = table

    def create_star_image(
        self,
        star_coordinates: tuple,
        ordinary_ray: float,
        extra_ordinary_ray: float = 0,
        seeing: float = 1,
        seed: float = 1,
    ) -> ndarray:
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

        seed: float, optional
            The seed used for the creation of the noise images. Default value to 1.

        Returns
        -------

        star_image: ndarray
            The Point Spred Function of the star for the respective CCD operation mode.
        """
        self.seed = seed
        self._create_table(star_coordinates, seeing)
        star_image = self._create_image_ordinary_ray(ordinary_ray)
        if extra_ordinary_ray != 0:
            star_image += self._create_image_extra_ordinary_ray(extra_ordinary_ray)
        return star_image

    def _create_image_ordinary_ray(self, ordinary_ray) -> ndarray:
        gaussian_amplitude = self._calculate_gaussian_amplitude(ordinary_ray)
        self.table["amplitude"] = gaussian_amplitude

        star_image = make_model_image(
            self.shape,
            Gaussian2D(),
            self.table,
            x_name="x_mean",
            y_name="y_mean",
        )
        star_image += self._make_noise_image()
        return star_image

    def _create_image_extra_ordinary_ray(self, extra_ordinary_ray) -> ndarray:
        gaussian_amplitude = self._calculate_gaussian_amplitude(extra_ordinary_ray)
        self.table["amplitude"] = gaussian_amplitude
        self.table["x_mean"] -= self.SPARC4_POL_SEPARATION
        self.table["y_mean"] -= self.SPARC4_POL_SEPARATION

        star_image = make_model_image(
            self.shape,
            Gaussian2D(),
            self.table,
            x_name="x_mean",
            y_name="y_mean",
        )
        star_image += self._make_noise_image()

        return star_image

    def _calculate_gaussian_amplitude(self, photons_per_second) -> float:
        t_exp = self.ccd_operation_mode["t_exp"]
        em_gain = self.ccd_operation_mode["em_gain"]
        binn = self.ccd_operation_mode["binn"]

        gaussian_amplitude = (
            photons_per_second
            * t_exp
            * em_gain
            * binn**2
            / (self.ccd_gain * 2 * pi * self.table["x_stddev"] * self.table["y_stddev"])
        )
        return gaussian_amplitude

    def _make_noise_image(self) -> ndarray:
        amplitude = self.table["amplitude"]
        noise_image = (
            make_noise_image(self.shape, "poisson", amplitude, seed=self.seed)
            - amplitude
        )
        self.table["amplitude"] = 1
        noise_image *= make_model_image(
            self.shape, Gaussian2D(), self.table, x_name="x_mean", y_name="y_mean"
        )
        return noise_image

    def calculate_npix_star(self, seeing: float) -> int:
        """Calculate the number of pixel of PSF of the object.

        Parameters
        ----------
            seeing (float): seeing of the object.

        Returns
        -------
            npix (int): number of pixels of PSF of the object.
        """
        gaussian_std = seeing / (
            self.SPARC4_PLATE_SCALE * self.ccd_operation_mode["binn"] * 2
        )
        fwhm = 2.355 * gaussian_std
        psf_radius = 3 * fwhm
        npix = pi * psf_radius**2
        return npix
