# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 13:26:15 2021

@author: denis
"""


from .AIS import Artificial_Image_Simulator


if __name__ == "__main__":

    ccd_operation_mode = dic = {'em_mode': 0, 'em_gain': 1, 'binn': 1,
                                'preamp': 1, 'hss': 1, 't_exp': 1,
                                'ccd_temp': -70}

    ais = Artificial_Image_Simulator(
        star_magnitude=100, sky_magnitude=10, gaussian_std=3,
        ccd_operation_mode=ccd_operation_mode, channel=1,
        image_dir=r'C:\Users\denis\Desktop\AIS')

    # ais.create_artificial_image()
