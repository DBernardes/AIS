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


@pytest.fixture
def source():
    return Source()


sed = np.linspace(100, 1000, 100)


def test_get_sed(source):
    source.write_sed(sed)
    assert np.allclose(source.get_sed(), sed)


magnitude = 10


def test_get_sed_error(source):
    with pytest.raises(ValueError):
        source.calculate_sed('error', magnitude)


# ------------------------------------------------------------
wv = np.linspace(350, 1100, 100)*1e-9
temperature = 5700
sed_blackbody = 2 * h * c ** 2 / wv ** 5 * 1 / \
    (np.exp(h * c / (wv * k * temperature)) - 1)


def test_calculate_sed_blackbody(source):
    assert np.allclose(source._calculate_sed_blackbody(
        wv, temperature), sed_blackbody)


TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
EFFECT_WAVELENGTH = 550e-9  # m
S_0 = vega_fluxd.get()['Johnson V'].value*1e7  # W/m2/m
effective_flux = S_0*10**(-magnitude/2.5) * \
    TELESCOPE_EFFECTIVE_AREA*EFFECT_WAVELENGTH/h*c


def test_calculate_effective_flux(source):
    assert source._calculate_effective_flux(magnitude) == effective_flux


# ------------------------------------------------------------

wavelength_interval = (350, 1100, 100)
calculation_method = 'blackbody'
spl = splrep(wv, sed_blackbody)
normalization_flux = splev(EFFECT_WAVELENGTH, spl)
new_sed = sed_blackbody * effective_flux / normalization_flux


def test_calculate_sed(source):
    class_sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature)
    assert np.allclose(class_sed, new_sed, rtol=0.005)
