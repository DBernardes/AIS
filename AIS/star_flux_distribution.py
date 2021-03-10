#!/usr/bin/env python
# coding: utf-8
#Author: Denis Varise Bernardes.
#Date: 04/04/2020. 

from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image
import numpy as np

class Star_Flux_Distribution:

    def __init__(self, ccd_info, star_flux, gaussian_stddev):        
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
        self.star_flux = star_flux
        self.gaussian_stddev = gaussian_stddev
        self.star_flux_distribution = []
    

    def create_star_flux_distribution(self):
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
        gaussian_stddev = self.gaussian_stddev

        #Calculation of the guassian amplitude based on the parameters provided to the software.        
        gaussian_amplitude = self.star_flux * t_exp * em_gain * binn**2 / gain
        #The create image has 200 x 200 pixels
        shape = (nlines, ncolumns)
        table = Table()
        table['amplitude'] = [gaussian_amplitude]
        #X coordinate of the star
        table['x_mean'] = [nlines//2]
        #Y coordinate of the star
        table['y_mean'] = [ncolumns//2]
        #Standard devion of the Gaussian in X direction
        table['x_stddev'] = [gaussian_stddev/binn]
        #Standard devion of the Gaussian in Y direction
        table['y_stddev'] = [gaussian_stddev/binn]
        #Rotation angle of the Gaussian. It is fixed in zero.
        table['theta'] = np.radians(np.array([0]))
        #This command creates the 2D-Gaussiain Distribution
        self.star_flux_distribution = make_gaussian_sources_image(shape, table)

        return self.star_flux_distribution
        


  
