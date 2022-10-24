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

from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength


@pytest.fixture
def source():
    return Source()


sed = np.linspace(100, 1000, 100)


calculation_method = 'blackbody'


def test_get_sed(source):
    source.write_sed(sed)
    assert np.allclose(source.get_sed(), sed)


wv = np.linspace(350, 1100, 100)*1e-9
temperature = 5700
magnitude = 10
effective_wavelength = 550
wavelength_interval = (350, 1100, 100)


def test_calculate_sed_spectral_standard(source):
    assert np.allclose(source.calculate_sed(
        'spectral_standard', wavelength_interval, magnitude, effective_wavelength), [])


def test_get_sed_error(source):
    with pytest.raises(ValueError):
        source.calculate_sed('error', wavelength_interval,
                             magnitude, effective_wavelength)

# ------------------------------------------------------------


sed_blackbody = 2 * h * c ** 2 / wv ** 5 * 1 / \
    (np.exp(h * c / (wv * k * temperature)) - 1)


def test_calculate_sed_blackbody(source):
    assert np.allclose(source._calculate_sed_blackbody(
        wv, temperature), sed_blackbody)


TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
S_0 = 4e-2  # W/m2/m
effective_flux = S_0*10**(-magnitude/2.5) * \
    TELESCOPE_EFFECTIVE_AREA*effective_wavelength*1e-9/h*c


def test_calculate_effective_flux(source):
    assert source._calculate_effective_flux(
        magnitude, effective_wavelength) == effective_flux


# ------------------------------------------------------------
new_sed = sed_blackbody * effective_flux / max(sed_blackbody)


def test_calculate_sed(source):
    assert np.allclose(source.calculate_sed(
        calculation_method, wavelength_interval, magnitude, effective_wavelength, temperature), new_sed)
