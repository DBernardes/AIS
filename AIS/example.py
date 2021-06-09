# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:23:27 2021

@author: denis
"""


from Artificial_Image_Simulator import Artificial_Image_Simulator

dic = {
    "em_mode": 0,
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 1024,
}


ais = Artificial_Image_Simulator(
    star_magnitude=15,
    sky_magnitude=10,
    ccd_operation_mode=dic,
    channel=1,
    gaussian_std=3,
    star_coordinates=[100, 100],
    bias_level=500,
    sparc4_operation_mode="phot",
    image_dir=r"C:\Users\denis\Desktop\FITS",
)

ais.create_artificial_image()
