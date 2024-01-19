"""

Background Image
================

This is the Background Image Class used to generate a background image.
"""


import os

import numpy as np
import pandas as pd
from numpy import ndarray
from photutils.datasets import make_noise_image

from ..Noise import Noise

__all__ = ["Background_Image"]


class Background_Image:
    """
    Background Image Class.

    Parameters
    ----------
    ccd_operation_mode: dictionary
        A python dictionary with the CCD operation mode.
        The allowed keywords values for the dictionary are

        em_mode : ['EM', 'Conv']
            Electron Multiplying mode of the camera.
        em_gain : float
            EM gain of the camera.
        readout : [0.1, 1, 10, 20, 30]
            Readout rate of the pixels in MHz.
        preamp : [1, 2]
            Pre-amplifier gain.
        binn : [1, 2]
            Binning of the pixels.
        t_exp : float
            Exposure time in seconds.
        image_size : int
            Image size in pixels.
        temp: float
            The camera temperature in celsius degree.
    channel: integer
        The channel related to the camera.
    ccd_temp: float
        The CCD temperature in celsius degree.
    bias_level : integer (optional)
        The bias level of the image in ADU.
    """

    PIXEL_SENSIBILITY = 0.03
    SPREADSHEET_PATH = os.path.join("AIS", "Background_Image", "preamp_gains.csv")

    def __init__(
        self,
        ccd_operation_mode: dict,
        channel: int,
        ccd_temp: float | int = -70,
        bias_level: int = 500,
    ) -> None:
        """Initialize the class."""
        self.bias_level = bias_level
        self.channel = channel
        self.ccd_operation_mode = ccd_operation_mode

        self._NOISE_FACTOR = 1
        if ccd_operation_mode["em_mode"] == "EM":
            self._NOISE_FACTOR = 1.4

        noise = Noise(channel)  # ? Injeção de dependência
        self.read_noise = noise.calculate_read_noise(ccd_operation_mode)
        self.dark_noise = (
            noise.calculate_dark_current(ccd_temp) * ccd_operation_mode["t_exp"]
        )
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

    def create_bias_background(self) -> ndarray:
        """Create the bias background.

        This functions creates a bias background with a noise distribution given by a gaussian distribution over the read noise.

        Returns
        -------
        bias_background: array like
            A bias background for the respective CCD operation mode.

        """

        image_size = self.ccd_operation_mode["image_size"]
        shape = (image_size, image_size)
        noise_adu = self.read_noise / self.ccd_gain
        bias_background = make_noise_image(
            shape, distribution="gaussian", mean=self.bias_level, stddev=noise_adu
        )

        return bias_background

    def create_dark_background(self) -> ndarray:
        """
        Create a dark background.

        This functions creates a dark background given by the ccd operation mode,
        the dc noise, and the bias level. This background is created woth a noise given by a gaussian
        distribution over the read noise and dc noise. Also, the
        extra noise of the EM amplification is considered.

        Returns
        -------
        dark_background: array like
            A dark background for the respective CCD operation mode and the dc level

        """

        image_size = self.ccd_operation_mode["image_size"]
        em_gain = self.ccd_operation_mode["em_gain"]
        binn = self.ccd_operation_mode["binn"]
        shape = (image_size, image_size)
        dark_level = (
            self.bias_level + self.dark_noise * em_gain * binn**2 / self.ccd_gain
        )

        noise = (
            np.sqrt(
                self.read_noise**2
                + self.dark_noise * self._NOISE_FACTOR**2 * em_gain**2 * binn**2
            )
            / self.ccd_gain
        )

        dark_background = make_noise_image(
            shape, distribution="gaussian", mean=dark_level, stddev=noise
        )

        return dark_background

    def create_flat_background(self) -> ndarray:
        """
        Create a flat background.

        This functions creates a flat background with a noise distribution given by the contribution of the
        the read noise, the Poisson noise, and the pixel sensibility noise.
        The extra noise related with the EM amplification is also considered.

        Returns
        -------
        flat_background: array like
            A flat background with counts distribution around half of the
            pixels depth.

        """
        em_gain = self.ccd_operation_mode["em_gain"]
        binn = self.ccd_operation_mode["binn"]
        image_size = self.ccd_operation_mode["image_size"]
        FLAT_LEVEL = 2**15  # ADU
        if self.ccd_operation_mode["preamp"] == 1:
            FLAT_LEVEL /= 2

        poisson_noise = FLAT_LEVEL / self.ccd_gain
        shape = (image_size, image_size)

        noise = (
            np.sqrt(
                self.read_noise**2
                + (poisson_noise * (1 + self.PIXEL_SENSIBILITY) + self.dark_noise)
                * self._NOISE_FACTOR**2
                * em_gain**2
                * binn**2
            )
            / self.ccd_gain
        )

        flat_background = make_noise_image(
            shape, distribution="gaussian", mean=FLAT_LEVEL, stddev=noise
        )

        return flat_background

    def create_sky_background(self, sky_flux: float) -> ndarray:
        """
        Create a sky background.

        This functions creates a sky background given
        by the sky flux, the dark noise, and the bias
        level. Over this background, there is a noise given by the read noise, dark noise, and sky noise.
        The extra noise of the EM amplification is also considered.

        Parameters
        ----------

        sky_flux : float
            Photons/s of the sky.

        Returns
        -------
        sky_background: array like
            A sky background for the respective sky flux, and the dark noise.

        """
        t_exp = self.ccd_operation_mode["t_exp"]
        em_gain = self.ccd_operation_mode["em_gain"]
        binn = self.ccd_operation_mode["binn"]
        image_size = self.ccd_operation_mode["image_size"]

        sky_background = (
            self.bias_level
            + (self.dark_noise + sky_flux) * t_exp * em_gain * binn**2 / self.ccd_gain
        )

        noise = (
            np.sqrt(
                self.read_noise**2
                + (sky_flux * t_exp + self.dark_noise)
                * self._NOISE_FACTOR**2
                * em_gain**2
                * binn**2
            )
            / self.ccd_gain
        )

        shape = (image_size, image_size)
        sky_background = make_noise_image(
            shape, distribution="gaussian", mean=sky_background, stddev=noise
        )

        return sky_background
