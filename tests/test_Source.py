# # -*- coding: utf-8 -*-
# Tests of the Source Class
# Oct 24th 2022
# @author: denis
#

import os
from copy import copy
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Energy_Distribution import Source
from scipy.interpolate import splev, splrep
from scipy.constants import c, h, k
from sbpy.calib import vega_fluxd
from math import pi, sqrt, tan
from sys import exit
from AIS.Spectral_Response._utils import calculate_polarizer_matrix, apply_matrix

obj_wavelength = np.linspace(400, 1100, 100)


@pytest.fixture
def source():
    return Source()


SPECTRAL_LIB_PATH = os.path.join(
    "AIS", "Spectral_Energy_Distribution", "Spectral_Library"
)


def test_spectral_lib_path(source):
    assert source.SPECTRAL_LIB_PATH == SPECTRAL_LIB_PATH


def test_get_sed(source):
    sed = np.linspace(100, 1000, 100)
    source.write_sed(sed)
    assert np.allclose(source.get_sed(), sed)


magnitude = 10


def test_get_sed_error(source):
    with pytest.raises(ValueError):
        source.calculate_sed("error", magnitude)


# ------------------------------------------------------------
wv = np.linspace(350, 1100, 100)
temperature = 5700
sed_blackbody = (
    2
    * pi
    * h
    * c**2
    / (wv * 1e-9) ** 5
    * 1
    / (np.exp(h * c / (wv * 1e-9 * k * temperature)) - 1)
)


def test_calculate_sed_blackbody(source):
    class_sed = source._calculate_sed_blackbody(wv, temperature)
    assert np.allclose(class_sed, sed_blackbody)


# ------------------------------------------------------------
TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
EFFECT_WAVELENGTH = 555.6  # nm
S_0 = 3.658e-2  # W/m2/m
effective_flux = (
    S_0
    * 10 ** (-magnitude / 2.5)
    * TELESCOPE_EFFECTIVE_AREA
    * EFFECT_WAVELENGTH
    * 1e-9
    / (h * c)
)
wavelength_interval = (350, 1100, 100)
calculation_method = "blackbody"
spl = splrep(wv, sed_blackbody)
normalization_flux = splev(EFFECT_WAVELENGTH, spl)
new_sed = sed_blackbody * effective_flux / normalization_flux


def test_calculate_sed(source):
    class_wv, class_sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature
    )
    temp = np.zeros((4, 100))
    temp[0] = new_sed
    assert np.allclose(class_wv, wv, rtol=0.005)
    assert np.allclose(temp, class_sed, rtol=0.005)


SPECTRAL_LIB_PATH = os.path.join(
    "AIS", "Spectral_Energy_Distribution", "Spectral_Library"
)
NAME_SED_SPECTRAL_TYPE = os.listdir(SPECTRAL_LIB_PATH)


def test_read_spectral_library(source):
    for name in NAME_SED_SPECTRAL_TYPE:
        name = name.split(".")[0][2:]
        if name == "":
            continue
        wv, sed = source._read_spectral_library(name)
        path = os.path.join(SPECTRAL_LIB_PATH, "uk" + name + ".csv")
        file_data = pd.read_csv(path)
        assert np.allclose(wv, file_data["wavelength (nm)"])
        assert np.allclose(sed, file_data["flux (F_lambda)"])


def test_get_calculate_sed_spectral_lib(source):
    for name in NAME_SED_SPECTRAL_TYPE:
        name = name.split(".")[0][2:]
        if name == "":
            continue
        wv, sed = source.calculate_sed(
            "spectral_library", magnitude, spectral_type=name
        )
        path = os.path.join(SPECTRAL_LIB_PATH, "uk" + name + ".csv")
        file_data = pd.read_csv(path)
        new_sed = file_data["flux (F_lambda)"] * effective_flux
        temp = np.zeros((4, len(new_sed)))
        temp[0] = new_sed
        assert np.allclose(wv, file_data["wavelength (nm)"], rtol=0.005)
        assert np.allclose(sed, temp, rtol=0.005)


def test_linear_polarization(source):
    percent_pol = 70
    pol_angle = 10

    _, sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature
    )
    polarized_sed = source.apply_linear_polarization(percent_pol, pol_angle)

    theta = np.deg2rad(pol_angle)
    tan_value = tan(2 * theta)
    sed[1] = sed[0] * percent_pol / (100 * sqrt(1 + tan_value**2))
    sed[2] = sed[1] * tan_value

    assert np.allclose(sed, polarized_sed, atol=1e-3)


def test_circular_polarization(source):
    percent_pol = 70

    _, sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature
    )
    polarized_sed = source.apply_circular_polarization(percent_pol)

    sed[3] = percent_pol * sed[0] / 100

    assert np.allclose(sed, polarized_sed, atol=1e-3)


def test_polarization(source):
    stokes = [0.1, 0.2, 0.3]
    _, sed = source.calculate_sed(
        calculation_method, magnitude, wavelength_interval, temperature
    )
    polarized_sed = source.apply_polarization(stokes)

    I = sed[0]
    sed[1] = I * stokes[0]
    sed[2] = I * stokes[1]
    sed[3] = I * stokes[2]

    assert np.allclose(sed, polarized_sed)
