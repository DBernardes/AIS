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
from scipy.constants import c, h, k
from sbpy.calib import vega_fluxd


from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength


@pytest.fixture
def sky():
    return Sky()


BASE_PATH = os.path.join('AIS', 'Spectral_Energy_Distribution')
csv_file = os.path.join(BASE_PATH, 'moon_magnitude.csv')
ss = pd.read_csv(csv_file)


def test_read_csv(sky):
    for value_name in ['new', 'waxing', 'waning', 'full']:
        wv, value = sky._read_csv('moon_magnitude.csv', value_name)
        assert np.allclose(wv, ss['wavelength'])
        assert np.allclose(value, ss[value_name])


moon_phase = 'new'
EFFECT_WAVELENGTH = 555.6  # nm
TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
S_0 = vega_fluxd.get()["Johnson V"].value*1e7  # W/(m.m2)
new_sed = S_0*10**(-ss[moon_phase]/2.5)*TELESCOPE_EFFECTIVE_AREA * \
    EFFECT_WAVELENGTH*1e-9/h*c
spl = splrep(ss['wavelength'], new_sed)
new_sed = splev(obj_wavelength, spl)


def test_calculate_sef(sky):
    sed = sky.calculate_sed(moon_phase, obj_wavelength)
    assert np.allclose(sed, new_sed)
