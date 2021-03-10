#!/usr/bin/env python
# coding: utf-8
#Author: Denis Varise Bernardes.
#Date: 04/04/2020.

import numpy as np
from photutils.datasets import make_noise_image

class Create_Noise_Image:

    def __init__(self, ccd_info, sky_flux):        
        self.t_exp = ccd_info['t_exp']
        self.em_gain = ccd_info['em_gain']
        self.gain = ccd_info['gain']
        self.bias_level = ccd_info['bias_level']
        self.dark_noise = ccd_info['dark_noise']
        self.read_noise = ccd_info['read_noise']
        self.noise_factor = 1
        if ccd_info['em_mode'] == 1:
            self.noise_factor = 1.4            
        self.bin = ccd_info['bin']
        self.lines_number = ccd_info['lines_number']
        self.columns_number = ccd_info['columns_number']
        self.sky_flux = sky_flux       
        self.noise_image = []


    def create_noise_image(self):
        #This function creates the artificial image.
        t_exp = self.t_exp
        em_gain = self.em_gain
        gain = self.gain
        bias = self.bias_level                
        dc = self.dark_noise
        rn = self.read_noise               
        nf = self.noise_factor
        binn = self.bin
        nlines = self.lines_number
        ncolumns = self.columns_number
        sky = self.sky_flux        

        shape = (nlines,ncolumns)        
        #Calculation of the background level of the image.        
        background_level = bias + (dc + sky) * t_exp * em_gain * binn**2 / gain
        #Calculation of the image noise
        image_noise = np.sqrt(rn**2 + (sky + dc)*t_exp * nf**2 * em_gain**2 * binn**2)/gain               
        #This command creates a noise image with a Gaussian Distribution. This image has counts distribution
        # arround zero, for the noise calculated in the previous step.
        self.noise_image = make_noise_image(shape, distribution='gaussian', mean=background_level, stddev=image_noise)

        return self.noise_image
        
        
        
        
