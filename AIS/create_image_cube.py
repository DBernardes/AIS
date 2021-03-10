#!/usr/bin/env python
# coding: utf-8
#14/10/2020. Denis Varise Bernardes.

import star_flux_distribution as SFD
import noise_image as NI
import create_header as CH
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import Read_Noise_Calc as RNC
import astropy.io.fits as fits

class Create_Image_Cube:

    def __init__(self, ccd_info, star_flux, sky_flux, gaussian_stddev = 3, image_dir='', image_name=''):
        if ccd_info['em_mode'] == 0: ccd_info['em_gain'] = 1

        try: ccd_info['bias_level']
        except: ccd_info['bias_level'] = 0

        try: ccd_info['lines_number']
        except: ccd_info['lines_number'] = 200

        try: ccd_info['columns_number']
        except: ccd_info['columns_number'] = 200

        try: ccd_info['serial_number']
        except: ccd_info['serial_number'] = 9914

        try: ccd_info['ccd_temperature']
        except: ccd_info['ccd_temperature'] = -70     
        
        self.ccd_info = ccd_info        
        self.star_flux = star_flux
        self.sky_flux = sky_flux
        self.image_dir = '..\\'
        if image_dir != '':        
            if '\\' not in image_dir[-1]: image_dir+='\\'
            self.image_dir = image_dir                
        #Name of the createde image. It is automatically generated.
        self.image_name = image_name
        self.write_image_name()

        #The software seeks for these three paramaters in the ccd_info dictionary.
        #If is does not find, it sets these parameters based on the CCD operation mode provided.
        self.set_dark_current()
        self.set_read_noise()
        self.set_gain()
        
        self.SFD = SFD.Star_Flux_Distribution(ccd_info, star_flux, gaussian_stddev)    
        self.NI = NI.Create_Noise_Image(ccd_info, sky_flux)
        self.CH = CH.Create_Image_Header(ccd_info, sky_flux, self.image_name)
        self.image_cube = np.zeros((ccd_info['cube_size'], ccd_info['lines_number'], ccd_info['columns_number']))


    def set_dark_current(self):
        #This function seeks for the dark current noise in the ccd_info dictionary.
        #If it does not find, it calculates the dark current of the CCD based on the model presentes in
        #https://ui.adsabs.harvard.edu/abs/2018PASP..130i5002B/abstract.
        serial_number = self.ccd_info['serial_number']
        T = self.ccd_info['ccd_temperature']
        try:self.ccd_info['dark_noise']
        except:
            dark_current = 0         
            if serial_number == 9914:
                dark_current = 24.66*np.exp(0.0015*T**2+0.29*T) 
            if serial_number == 9915:
                dark_current = 35.26*np.exp(0.0019*T**2+0.31*T)
            if serial_number == 9916:
                dark_current = 9.67*np.exp(0.0012*T**2+0.25*T)
            if serial_number == 9917:
                dark_current = 5.92*np.exp(0.0005*T**2+0.18*T)
            self.ccd_info['dark_noise'] = dark_current*self.ccd_info['t_exp']
        
        

    def set_read_noise(self):
        #If the read noise was not provided to the ccd_info parameter,
        #this function calls the ReadNoiseCalc library to calculate the read noise
        # of the CCD, based on the provided operation mode
        try:self.ccd_info['read_noise']
        except:        
            RN = RNC.ReadNoiseCalc()    
            RN.write_operation_mode(self.ccd_info['em_mode'],
                                    self.ccd_info['em_gain'],
                                    self.ccd_info['hss'],
                                    self.ccd_info['preamp'],
                                    self.ccd_info['bin'])
            RN.calc_read_noise()        
            self.ccd_info['read_noise'] = RN.calc_read_noise()
        


    def set_gain(self):
        #This function seeks for the gain value in the ccd_info dictionary provided.
        #If it does not find, it  sets the CCD gain based on the provided operation mode.
        #The gain values in this function were obtained through the SPARC4's camera datasheet
        try: self.ccd_info['gain']
        except:        
            em_mode = self.ccd_info['em_mode']
            hss = self.ccd_info['hss']
            preamp = self.ccd_info['preamp']
            gain = 0
            if em_mode == 1:
                if hss == 30:
                    if preamp == 1:
                        gain = 17.2
                    if preamp == 2:
                        gain = 5.27
                if hss == 20:
                    if preamp == 1:
                        gain = 16.4
                    if preamp == 2:
                        gain = 4.39
                if hss == 10:
                    if preamp == 1:
                        gain = 16.0
                    if preamp == 2:
                        gain = 3.96
                if hss == 1:
                    if preamp == 1:
                        gain = 15.9
                    if preamp == 2:
                        gain = 3.88
            else:
                if hss == 1:
                    if preamp == 1:
                        gain = 3.37
                    if preamp == 2:
                        gain = 0.8
                if hss == 0.1:
                    if preamp == 1:
                        gain = 3.35
                    if preamp == 2:
                        gain = 0.8
            self.ccd_info['gain'] = gain        




    def write_image_name(self):
        #If the image name is not provided,
        #this function creates the image name based on the provided operation mode
        if self.image_name == '':
            em_gain = '_G' + str(self.ccd_info['em_gain'])
            em_mode = 'CONV_'
            if self.ccd_info['em_mode'] == 1: em_mode = 'EM_'
            hss = str(self.ccd_info['hss']) + 'MHz'
            preamp = '_PA' + str(self.ccd_info['preamp'])
            binn = '_B' + str(self.ccd_info['bin'])
            t_exp = '_TEXP'+ str(self.ccd_info['t_exp'])             
            self.image_name = em_mode + hss + preamp + binn + t_exp + em_gain
        

        


    def create_image_cube(self):
        star_image = self.SFD.create_star_flux_distribution()
        noise_image = self.NI.create_noise_image()        
        for i in range(self.ccd_info['cube_size']):
            self.image_cube[i] = star_image + noise_image


    def save_image_cube(self):
        hdr = self.CH.create_image_header()
        fits.writeto(self.image_dir + self.image_name+ '.fits',
                     self.image_cube,
                     overwrite=True,
                     header=hdr)       
        
   
