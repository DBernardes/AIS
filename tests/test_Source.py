# # -*- coding: utf-8 -*-
# Tests of the Source Class
# Oct 24th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Energy_Distribution import Source
from scipy.interpolate import splev, splrep
from scipy.constants import c, h, k
from sbpy.calib import vega_fluxd
from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength
from math import pi
from sys import exit


@pytest.fixture
def source():
    return Source()


SPECTRAL_LIB_PATH = os.path.join(
    'AIS', 'Spectral_Energy_Distribution', 'Spectral_Library')


def test_spectral_lib_path(source):
    assert source.SPECTRAL_LIB_PATH == SPECTRAL_LIB_PATH


def test_get_sed(source):
    sed = np.linspace(100, 1000, 100)
    source.write_sed(sed)
    assert np.allclose(source.get_sed(), sed)


magnitude = 10


def test_get_sed_error(source):
    with pytest.raises(ValueError):
        source.calculate_sed('error', magnitude)


# ------------------------------------------------------------
wv = np.linspace(350, 1100, 100)
temperature = 5700
sed_blackbody = 2 * pi * h * c ** 2 / (wv * 1e-9) ** 5 * 1 / \
    (np.exp(h * c / (wv * 1e-9 * k * temperature)) - 1)


def test_calculate_sed_blackbody(source):
    class_sed = source._calculate_sed_blackbody(
        wv, temperature)
    assert np.allclose(class_sed, sed_blackbody)


# ------------------------------------------------------------
TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
EFFECT_WAVELENGTH = 550  # nm
S_0 = vega_fluxd.get()['Johnson V'].value*1e7  # W/m2/m
effective_flux = S_0*10**(-magnitude/2.5) * \
    TELESCOPE_EFFECTIVE_AREA*EFFECT_WAVELENGTH*1e-9/h*c
wavelength_interval = (350, 1100, 100)
calculation_method = 'blackbody'
spl = splrep(wv, sed_blackbody)
normalization_flux = splev(EFFECT_WAVELENGTH, spl)
new_sed = sed_blackbody * effective_flux / normalization_flux


def test_calculate_sed(source):
    class_wv, class_sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature)
    assert np.allclose(class_wv, wv)
    assert np.allclose(class_sed, new_sed)


SPECTRAL_LIB_PATH = os.path.join(
    'AIS', 'Spectral_Energy_Distribution', 'Spectral_Library')
NAME_SED_SPECTRAL_TYPE = os.listdir(SPECTRAL_LIB_PATH)


def test_read_spectral_library(source):
    for name in NAME_SED_SPECTRAL_TYPE:
        name = name.split('.')[0][2:]
        if name == '':
            continue
        wv, sed = source._read_spectral_library(name)
        path = os.path.join(SPECTRAL_LIB_PATH, 'uk' + name + '.dat')
        file_data = np.loadtxt(path)
        assert np.allclose(wv, file_data[:, 0]/10)
        assert np.allclose(sed, file_data[:, 1])


def test_get_calculate_sed_spectral_lib(source):
    for key, val in NAME_SED_SPECTRAL_TYPE.items():
        wv, sed = source.calculate_sed(
            'spectral_library', magnitude, spectral_type=key)
        path = os.path.join(SPECTRAL_LIB_PATH, val)
        file_data = np.loadtxt(path)
        assert np.allclose(wv, file_data[:, 0]/10, rtol=0.005)
        assert np.allclose(sed, file_data[:, 1]*effective_flux, rtol=0.005)
