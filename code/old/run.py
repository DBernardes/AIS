#!/usr/bin/env python
# coding: utf-8
"""
This file runs the Artificial Images Generator algorithm
For more information, access https://github.com/DBernardes/Artificial-Images-Generator
"""
# 14/10/2020. Denis Varise Bernardes.
import matplotlib.pyplot as plt
import create_image_cube as CIC

ccd_info = {'t_exp': 1,
            'hss': 1,
            'em_mode': 1,
            'em_gain': 2,
            'bin': 1,
            'preamp': 2,
            'cube_size': 1}

star_flux = 2000
sky_flux = 12.39

CIC = CIC.Create_Image_Cube(ccd_info, star_flux, sky_flux)
CIC.create_image_cube()
CIC.save_image_cube()
