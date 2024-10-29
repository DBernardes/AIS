# # -*- coding: utf-8 -*-
# Tests of the Source Class
# Oct 24th 2022
# @author: denis
#

import os
import unittest
from copy import copy
from math import pi, sqrt, tan
from sys import exit

import numpy as np
import pandas as pd
import pytest
from sbpy.calib import vega_fluxd
from scipy.constants import c, h, k
from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit

from AIS.Spectral_Energy_Distribution import Source
from AIS.Spectral_Response._utils import apply_matrix, calculate_polarizer_matrix


class Test_Source(unittest.TestCase):
    EFFECT_WAVELENGTH = 545  # nm
    BASE_PATH = os.path.join("AIS", "Spectral_Energy_Distribution", "Spectral_Library")
    CSV_FILE = "moon_magnitude.csv"
    SOURCE = Source()
    MAGNITUDE = 10

    def setUp(self):
        self.WAVELENGTH = np.linspace(350, 1100, 100)
        return

    def test_calculate_sed_blackbody(self):
        wavelength = self.WAVELENGTH * 1e-9
        numerator = 2 * h * c**2
        denominator = wavelength**5 * (np.exp(h * c / (wavelength * k * 5700)) - 1)
        sed = numerator / denominator
        new_sed = self.SOURCE._calculate_sed_blackbody(self.WAVELENGTH, 5700)
        assert np.allclose(sed, new_sed)

    def test_read_spectral_library(self):
        _type = "a0i"
        path = os.path.join(self.BASE_PATH, "uk" + _type + ".csv")
        ss = pd.read_csv(path)
        ss["wavelength (nm)"], ss["flux (F_lambda)"]

        wv, sed = self.SOURCE._read_spectral_library(_type)

        assert np.allclose(wv, ss["wavelength (nm)"])
        assert np.allclose(sed, ss["flux (F_lambda)"])

    def test_read_spectral_library_error(self):
        with pytest.raises(FileNotFoundError):
            self.SOURCE._read_spectral_library("AAA")

    def test_print_available_spectral_types(self):
        self.SOURCE.print_available_spectral_types()

    def test_calculate_sed_bb(self):
        wv, sed = self.SOURCE.calculate_sed(
            "blackbody", self.MAGNITUDE, (350, 1100, 100), 5700
        )
        tmp = self.SOURCE._calculate_sed_blackbody(self.WAVELENGTH, 5700)
        normalization_flux = self.SOURCE._interpolate_spectral_distribution(
            self.WAVELENGTH, tmp, self.EFFECT_WAVELENGTH
        )
        tmp /= normalization_flux
        photons_density = self.SOURCE._calculate_photons_density(self.MAGNITUDE)

        new_sed = np.zeros((4, tmp.shape[0]))
        new_sed[0] = tmp * photons_density

        assert np.allclose(sed, new_sed)
        assert np.allclose(wv, self.WAVELENGTH)

    def test_calculate_sed_spec_lib(self):
        wv, sed = self.SOURCE.calculate_sed(
            "spectral_library", self.MAGNITUDE, (400, 1100, 100), spectral_type="a0i"
        )
        new_wv, tmp = self.SOURCE._read_spectral_library("a0i")
        tmp = self.SOURCE._interpolate_spectral_distribution(new_wv, tmp, wv)
        photons_density = self.SOURCE._calculate_photons_density(self.MAGNITUDE)

        new_sed = np.zeros((4, tmp.shape[0]))
        new_sed[0] = tmp * photons_density

        assert np.allclose(sed, new_sed)

    def test_linear_polarization(self):
        percent_pol = 70
        pol_angle = 10

        _, sed = self.SOURCE.calculate_sed(
            "blackbody", self.MAGNITUDE, (350, 1100, 100), 5700
        )

        theta = np.deg2rad(pol_angle)
        tan_value = tan(2 * theta)
        sed[1] = sed[0] * percent_pol / (100 * sqrt(1 + tan_value**2))
        sed[2] = sed[1] * tan_value

        polarized_sed = self.SOURCE.apply_linear_polarization(percent_pol, pol_angle)
        assert np.allclose(sed, polarized_sed, atol=1e-3)

    def test_circular_polarization(self):
        percent_pol = 70
        _, sed = self.SOURCE.calculate_sed(
            "blackbody", self.MAGNITUDE, (400, 1100, 100), 5700
        )
        sed[3] = percent_pol * sed[0] / 100
        polarized_sed = self.SOURCE.apply_circular_polarization(percent_pol)

        assert np.allclose(sed, polarized_sed, atol=1e-3)

    def test_polarization(self):
        stokes = [0.1, 0.2, 0.3]
        _, sed = self.SOURCE.calculate_sed(
            "blackbody", self.MAGNITUDE, (400, 1100, 100), 5700
        )

        I = sed[0]
        sed[1] = I * stokes[0]
        sed[2] = I * stokes[1]
        sed[3] = I * stokes[2]

        polarized_sed = self.SOURCE.apply_polarization(stokes)
        assert np.allclose(sed, polarized_sed)

    def test_Serkowski_curve(self):
        wavelength, p_max, K, l_max = 450e-9, 10, 1.15, 550e-9
        p = p_max * np.exp(-K * np.log(l_max / wavelength) ** 2)
        assert self.SOURCE._Serkowski_curve(wavelength, p_max, l_max) == p

    def test_verify_pol_vals(self):
        with pytest.raises(ValueError):
            self.SOURCE._verify_pol_BVRI({"B": 1, "V": 1, "R": 1, "I": None})
        with pytest.raises(ValueError):
            self.SOURCE._verify_pol_BVRI({"B": 1, "V": 1, "R": 1, "I": -1})

    def test_adjust_Serkowski_curve(self):
        pol_BVRI = {"B": 7.812, "V": 7.029, "R": 5.951, "I": 4.706}
        popt, _ = curve_fit(
            self.SOURCE._Serkowski_curve,
            list(self.SOURCE.effect_wl.values()),
            list(pol_BVRI.values()),
            p0=(10, 500e-9),
        )
        assert np.allclose(popt, self.SOURCE._adjust_Serkowski_curve(pol_BVRI))

    def test_apply_Serkowski_curve(self):
        pol_BVRI = {"B": 7.812, "V": 7.029, "R": 5.951, "I": 4.706}
        sed = np.ones((4, 100))
        wv = np.linspace(400, 1100, 100)
        self.SOURCE.sed = sed
        self.SOURCE.wavelength = wv

        popt = self.SOURCE._adjust_Serkowski_curve(pol_BVRI)
        q_Stokes = self.SOURCE._Serkowski_curve(wv * 1e-9, *popt) / 100
        sed[1] = sed[0] * q_Stokes

        assert np.allclose(self.SOURCE.apply_Serkowski_curve(pol_BVRI), sed)
