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
from scipy.interpolate import splev, splrep, interp1d
from scipy.constants import c, h
from sbpy.calib import vega_fluxd


obj_wavelength = np.linspace(400, 1100, 100)


@pytest.fixture
def sky():
    return Sky()


BASE_PATH = os.path.join('AIS', 'Spectral_Energy_Distribution')
csv_file = os.path.join(BASE_PATH, 'moon_magnitude.csv')
ss = pd.read_csv(csv_file)


def test_read_csv(sky):
    for value_name in ['new', 'first quarter', 'third quarter', 'full']:
        wv, value = sky._read_csv('moon_magnitude.csv', value_name)
        assert np.allclose(wv, ss['wavelength'])
        assert np.allclose(value, ss[value_name])


moon_phase = 'new'
EFFECT_WAVELENGTH = 555.6  # nm
TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
S_0 = 3.658e-2  # W/(m.m2)
new_sed = S_0*10**(-ss[moon_phase]/2.5)*TELESCOPE_EFFECTIVE_AREA * \
    EFFECT_WAVELENGTH*1e-9/(h*c)
spl = splrep(ss['wavelength'], new_sed)
new_sed = np.zeros((4, 100))
new_sed[0] = splev(obj_wavelength, spl)


def test_calculate_sed(sky):
    sed = sky.calculate_sed(moon_phase, obj_wavelength)
    assert np.allclose(sed, new_sed)
