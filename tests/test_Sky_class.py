# # -*- coding: utf-8 -*-
# Tests of the Sky Class
# Oct 24th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Energy_Distribution import Sky
from scipy.interpolate import splev, splrep


from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength


@pytest.fixture
def sky():
    return Sky()


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
