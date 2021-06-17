# -*- coding: utf-8 -*-
"""Test of the Telescope Spectral Response Class.

Jun 14, 2021

@author: denis
"""

import numpy as np
import pytest
from AIS.Telescope_Spectral_Response import Telescope_Spectral_Response

l_init, l_final, l_step = 350, 1150, 50
wavelength_interv = np.asarray(range(l_init, l_final, l_step))
n = len(wavelength_interv)
specific_flux = np.ones((4, n))


@pytest.fixture
def tel_sr():
    return Telescope_Spectral_Response()


def test_apply_telescope_spectral_response(tel_sr):
    new_specific_flux = tel_sr.apply_telescope_spectral_response(
        specific_flux, l_init, l_final, l_step
    )
    assert np.allclose(specific_flux, new_specific_flux)
