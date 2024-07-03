"""Artificial Image Simulator Package.

This package was developed to create cubes of artificial star images,
simulating those images that would be acquired by using the SPARC4 1
CCD cameras in astronomical observations. To create the images, the AIS
models the star flux distribution as a 2D-Gaussian Distribution. This result is
added to a noise image, created based on the noise information of the SPARC4
CCDs, as a function of their operation mode.
"""

from .artificial_image_simulator import Artificial_Image_Simulator

__all__ = ["artificial_image_simulator"]
