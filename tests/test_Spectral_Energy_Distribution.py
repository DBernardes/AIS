# # -*- coding: utf-8 -*-
# Tests of the Spectral Energy Distributiion Class
# Oct 24th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Energy_Distribution import Spectral_Energy_Distribution
from scipy.interpolate import splev, splrep
from sbpy.calib import vega_fluxd
from scipy.constants import c, h, k

from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength


@pytest.fixture
def sed_obj():
    return Spectral_Energy_Distribution()


sed = np.linspace(100, 1000, 100)
wv = np.linspace(350, 1100, 100)


def test_get_sed(sed_obj):
    assert np.allclose(sed_obj.get_sed(), sed)


def test_write_sed(sed_obj):
    sed_obj.write_sed(sed)
    assert np.allclose(sed_obj.sed, sed)


def test_interpolate(sed_obj):
    spl = splrep(wv, sed)
    interpolated_sed = splev(obj_wavelength, spl)
    class_interpolated_sed = sed_obj._interpolate_spectral_distribution(
        wv, sed, obj_wavelength)
    assert np.allclose(class_interpolated_sed, interpolated_sed)


def test_calc_sed(sed_obj):
    assert np.allclose(sed_obj.calculate_sed(), sed)


magnitude = 10
TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
EFFECT_WAVELENGTH = 550  # nm
S_0 = vega_fluxd.get()['Johnson V'].value*1e7  # W/m2/m
effective_flux = S_0*10**(-magnitude/2.5) * \
    TELESCOPE_EFFECTIVE_AREA*EFFECT_WAVELENGTH*1e-9/h*c


def test_calculate_effective_flux(sed_obj):
    assert sed_obj._calculate_effective_flux(magnitude) == effective_flux
