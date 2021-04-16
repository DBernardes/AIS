"""
Artificial Images Simulator Package
===================================

The Artificial Images Simulator (AIS) package was developed to generate
artificial star images, similar to those images that would be acquired by
using the acquisition system of the instrument. To accomplish this,
the AIS models as star flux as a 2D gaussian distribution. Then, the star
flux is added to an image with a background level given by counts distribution
of an image of the SPARC4 cameras, as a function of its operation mode.
"""

# import read_noise_calc as RNC
# from photutils.datasets import make_noise_image
# import astropy.io.fits as fits
# import numpy as np
# import openpyxl

# from astropy.table import Table
# from photutils.datasets import make_gaussian_sources_image
# from sys import exit


__all__ = ['Artificial_Images_Simulator']


class Artificial_Images_Simulator:
    pass
