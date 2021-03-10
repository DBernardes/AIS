#!/usr/bin/env python
# coding: utf-8
#Author: Denis Varise Bernardes.
#Date: 04/04/2020.


import astropy.io.fits as fits
from sys import exit
import datetime


class Create_Image_Header:

    def __init__(self, ccd_info, sky_flux, image_name):        
        self.t_exp = ccd_info['t_exp']
        self.hss = ccd_info['hss']
        self.em_gain = ccd_info['em_gain']
        self.gain = ccd_info['gain']
        self.bias_level = ccd_info['bias_level']
        self.dark_noise = ccd_info['dark_noise']
        self.read_noise = ccd_info['read_noise']        
        self.em_mode = ccd_info['em_mode']
        self.bin = ccd_info['bin']
        self.lines_number = ccd_info['lines_number']
        self.columns_number = ccd_info['columns_number']
        self.image_name = image_name
        self.serial_number = ccd_info['serial_number']
        self.ccd_temperature = ccd_info['ccd_temperature']
        self.cube_size = ccd_info['cube_size']
        try: self.preamp = ccd_info['preamp']
        except: self.preamp = 2
        self.sky_flux = sky_flux        
        self.hdr = {}




    def create_image_header(self):               
        #This function creates the image header based on the paramters provided to the class.
        now = datetime.datetime.now()                
        date = str(now).split(' ')        
        date = date[0]+'T'+date[1]         
        
        hdr = fits.Header()                               
        hdr['NAXIS1']  =(self.lines_number, 'length of data axis 1')
        hdr['NAXIS2']  =(self.columns_number, 'length of data axis 2')                          
        hdr['EXTEND']  = ('T', 'FITS dataset may contain extensions')                                             
        hdr['COMMENT'] = 'and Astrophysics, volume 376, page 359'
        hdr['ACQMODE'] = ('Kinetic  ', 'Acquisition Mode')
        hdr['READMODE']= ('Image   ', 'Readout Mode')
        hdr['IMGRECT'] = ('1, '+ str(self.lines_number) + ' ' + str(self.columns_number) + ', 1', 'Image Format')
        hdr['HBIN']    = (self.bin, 'Horizontal Binning')
        hdr['VBIN']    = (self.bin,'Vertical Binning')
        hdr['INITLIN']    = ( 1 ,'Initial pixels line')
        hdr['INITCOL']    = ( 1 ,'Initial pixels column')
        hdr['FINALLIN']    = (self.lines_number,'Final pixels line')
        hdr['FINALCOL']    = (self.columns_number ,'Final pixels column')
        hdr['TRIGGER'] = ('External Start','Trigger Mode')
        hdr['EXPOSURE']= (self.t_exp, 'Total Exposure Time')
        hdr['TEMP']    = (self.ccd_temperature, 'Temperature')
        hdr['READTIME']= (str(1/self.hss)+'E-006' ,'Pixel readout time ')
        hdr['VSHIFT']  = ('0.6E-06', 'Vertical Shift Speed')
        hdr['GAIN']    = (self.gain, 'Preamp Gain (e-/ADU)')
        em_mode = 'Conventional'
        if self.em_mode == 1: em_mode = 'Electron Multiplying'
        hdr['OUTPTAMP']= (em_mode, 'Output Amplifier')
        hdr['EMGAIN'] = (self.em_gain, 'Electron Multiplying Gain')
        hdr['PREAMP']  = (str(self.preamp)+'x', 'Pre Amplifier Gain')
        hdr['SERNO']   = (self.serial_number, 'Serial Number')
        hdr['DATE']   = (date, 'File Creation Date (YYYY-MM-HHThh:mm:ss)')        
        hdr['IMAGE']  = (self.image_name, 'Nome do arquivo')
        hdr['IMG_DL0'] = ("{:.6e}".format(0), 'Delay da imagem 0')
        delay = 0
        n_pixels = self.lines_number*self.columns_number
        for i in range(self.cube_size-1):
            i+=1
            delay += (self.t_exp + n_pixels/(self.hss*1e6))
            hdr['IMG_DL'+str(i)] = ("{:.6e}".format(delay), 'Delay da imagem ' + str(i))
        self.hdr = hdr

        return self.hdr
