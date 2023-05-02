# -*- coding: utf-8 -*-
"""
Test of the Atmosphere Spectral Response Class.

Jun 14, 2021

@author: denis
"""


from AIS.Spectral_Response import Atmosphere
import pytest
import os
import numpy as np
from copy import copy
from scipy.optimize import curve_fit
from scipy.interpolate import splev, splrep, interp1d
import pandas as pd


obj_wavelength = np.linspace(400, 1100, 100)
_BASE_PATH = os.path.join('AIS', 'Spectral_Response')
csv_file_name = 'atmosphere_profile.csv'
_ATM_EXTINCTION = 'atmosphere.csv'
air_mass = 1
sky_condition = 'photometric'

atm_profile_ss = os.path.join(_BASE_PATH, csv_file_name) 
ss = pd.read_csv(atm_profile_ss)
wavelength, spectral_response = ss['Wavelength (nm)'], ss['Transmitance (%)']/100

def _func(x, c):
    csv = os.path.join(_BASE_PATH, csv_file_name)
    tmp_ss = pd.read_csv(csv)
    sys_wavelength, spectral_response = tmp_ss['Wavelength (nm)'], tmp_ss['Transmitance (%)']/100
    spl = interp1d(sys_wavelength, spectral_response,
                    bounds_error=False, kind='cubic')
    return spl(x)*c




@pytest.fixture
def atm():
    return Atmosphere()


def test_csv_file_name(atm):
    assert atm._CSV_FILE_NAME == csv_file_name


def test_base_path(atm):
    assert atm._BASE_PATH == _BASE_PATH

# --------------------------------------------------------------------------------------------


def test_read_csv_file(atm):
    new_wavelength, new_transmitance = atm._read_csv_file(
        atm_profile_ss)
    assert np.allclose(new_wavelength, wavelength)
    assert np.allclose(new_transmitance, spectral_response)

spl = interp1d(wavelength, spectral_response,
                       bounds_error=False, fill_value=0, kind='cubic')

def test_interpolate_spectral_response(atm):
    atm.obj_wavelength = obj_wavelength
    class_spectral_response = atm._interpolate_spectral_response(
        wavelength, spectral_response)    
    assert np.allclose(class_spectral_response, spl(obj_wavelength))

extinction_ss = pd.read_csv(os.path.join(_BASE_PATH, _ATM_EXTINCTION))        
extinction_coef = 10**(-0.4 * air_mass * extinction_ss[sky_condition])
popt_M1, _ = curve_fit(_func, extinction_ss['Wavelength (nm)'], extinction_coef) 

def test_get_spectral_response(atm):
    class_spectral_response = atm.get_spectral_response(
        obj_wavelength, air_mass, sky_condition)       
    new_spectral_response = spl(obj_wavelength)*popt_M1[0]        
        
    assert np.allclose(class_spectral_response, new_spectral_response)





def test_apply_spectral_response(atm):
    sed = np.ones((4, len(obj_wavelength)))
    class_reduced_esd = atm.apply_spectral_response(
        obj_wavelength, copy(sed), air_mass, sky_condition)
    
    new_spectral_response = spl(obj_wavelength)*popt_M1[0]  
    new_sed = np.multiply(new_spectral_response, sed)
    
    assert np.allclose(class_reduced_esd, new_sed)
