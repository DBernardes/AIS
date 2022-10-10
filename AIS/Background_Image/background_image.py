"""

Background Image
================

This is the Background Image Class used to generate a background image.
"""


import numpy as np
from photutils.datasets import make_noise_image


class Background_Image:
    """
    Background Image Class.

    Parameters
    ----------
    ccd_operation_mode: dictionary
        A python dictionary with the CCD operation mode.
        The allowed keywords values for the dictionary are

        em_mode : [0, 1]
            Electron Multiplying mode of the CCD.
        em_gain : float
            CCD Electron Multiplying gain
        hss : [0.1, 1, 10, 20, 30]
            Horizontal Shift Speed of the pixels
        preamp : [1, 2]
            Pre-amplifier gain
        binn : [1, 2]
            Binning of the pixels
        t_exp : float
            Exposure time in seconds
        image_size : int
            Image size in pixels

    ccd_gain : float
        CCD gain in e-/ADU.
    dark_current : float
        Dark current in e-/s of the CCD
    read_noise :
        Read noise in electrons of the CCD
    bias_level : integer
        The bias level of the image in ADU
    """

    _PIXEL_SENSIBILITY = 0.03

    def __init__(
        self,
        ccd_operation_mode,
        ccd_gain,
        dark_current,
        read_noise,
        bias_level,
    ):
        """Initialize the Background Image class."""
        self.ccd_gain = ccd_gain
        self.dark_current = dark_current
        self.read_noise = read_noise
        self.bias_level = bias_level

        self.em_mode = ccd_operation_mode["em_mode"]
        self.em_gain = ccd_operation_mode["em_gain"]
        self.binn = ccd_operation_mode["binn"]
        self.t_exp = ccd_operation_mode["t_exp"]
        self.preamp = ccd_operation_mode["preamp"]
        self.hss = ccd_operation_mode["hss"]
        self.image_size = ccd_operation_mode["image_size"]

        self._NOISE_FACTOR = 1.0
        if self.em_mode == 1:
            self._NOISE_FACTOR = 1.4

    def create_background_image(self, sky_flux):
        """
        Create a background image.

        This functions creates a background image with a background level given
        by the ccd operation mode, the sky flux, the dc noise, and the bias
        level. Over this image there is a noise given by a gaussian
        distribution over the read noise, dc noise, and sky noise. Also, the
        extra noise of the EM amplification is considered.

        Parameters
        ----------

        sky_flux : float
            Flux in photons/s of the sky

        Returns
        -------
        background_image: array like
            A background image for the respective CCD operation mode, the sky
            flux, and the dark current noise.

        """
        t_exp = self.t_exp
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        bias = self.bias_level
        dc = self.dark_current * t_exp
        rn = self.read_noise
        sky = sky_flux
        nf = self._NOISE_FACTOR
        binn = self.binn
        image_size = self.image_size

        shape = (image_size, image_size)
        background_level = bias + (dc + sky) * \
            t_exp * em_gain * binn ** 2 / ccd_gain

        noise = (
            np.sqrt(rn ** 2 + (sky + dc) * t_exp *
                    nf ** 2 * em_gain ** 2 * binn ** 2)
            / ccd_gain
        )
        self.background_image = make_noise_image(
            shape, distribution="gaussian", mean=background_level, stddev=noise
        )

        return self.background_image

    def create_dark_image(self):
        """
        Create a dark image.

        This functions creates a dark image with a background level given
        by the ccd operation mode, the dc noise, and the bias
        level. Over this image there is a noise given by a gaussian
        distribution over the read noise and dc noise. Also, the
        extra noise of the EM amplification is considered.

        Returns
        -------
        dark_image: array like
            A dark image for the respective CCD operation mode and the dc level

        """
        t_exp = self.t_exp
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        bias = self.bias_level
        dc = self.dark_current * t_exp
        rn = self.read_noise
        nf = self._NOISE_FACTOR
        binn = self.binn
        image_size = self.image_size

        shape = (image_size, image_size)
        dark_level = bias + (dc) * t_exp * em_gain * binn ** 2 / ccd_gain

        noise = (
            np.sqrt(rn ** 2 + (dc) * t_exp * nf **
                    2 * em_gain ** 2 * binn ** 2)
            / ccd_gain
        )

        self.dark_image = make_noise_image(
            shape, distribution="gaussian", mean=dark_level, stddev=noise
        )

        return self.dark_image

    def create_bias_image(self):
        """
        Create the bias image.

        This functions creates a bias image with the provided bias level given
        as a function of the ccd operation mode. Over this image there is a
        noise given by a gaussian distribution over the read noise.


        Returns
        -------
        bias_image: array like
            A bias image for the respective CCD operation mode.

        """
        ccd_gain = self.ccd_gain
        bias = self.bias_level
        rn = self.read_noise
        image_size = self.image_size

        shape = (image_size, image_size)
        noise = rn / ccd_gain
        self.bias_image = make_noise_image(
            shape, distribution="gaussian", mean=bias, stddev=noise
        )

        return self.bias_image

    def create_flat_image(self):
        """
        Create a flat image.

        This functions creates a flat image with a background level of
        half of the CCD's pixel depth. Over this background, the read noise,
        the Poisson noise, and the pixel sensibility noise are summed.
        Also, the extra noise of the EM amplification is considered.

        Returns
        -------
        flat_image: array like
            A flat image with counts distribution around half of the
            pixels depth.

        """
        em_gain = self.em_gain
        ccd_gain = self.ccd_gain
        BIAS = 32000  # ADU
        rn = self.read_noise
        nf = self._NOISE_FACTOR
        poisson_noise = BIAS / ccd_gain
        pixel_sensibility_noise = BIAS / ccd_gain * self._PIXEL_SENSIBILITY
        binn = self.binn
        image_size = self.image_size

        shape = (image_size, image_size)
        background_level = BIAS
        noise = (
            np.sqrt(
                rn ** 2
                + (poisson_noise + pixel_sensibility_noise)
                * nf ** 2
                * em_gain ** 2
                * binn ** 2
            )
            / ccd_gain
        )

        self.flat_image = make_noise_image(
            shape, distribution="gaussian", mean=background_level, stddev=noise
        )

        return self.flat_image
