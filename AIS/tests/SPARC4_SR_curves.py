# This file has all the spectral response curves used in the test_SR_SR.py file

import os
from sys import exit

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep

wavelength_interval = range(350, 1150, 50)
n = len(wavelength_interval)
specific_flux = np.ones((4, n))


def calculate_spline(transmitance, component_wavelength_interv, wavelength_interval):
    spl = splrep(component_wavelength_interv, transmitance)
    transmitance = splev(wavelength_interval, spl)
    return transmitance


# ------------------------------------ CCD ------------------------------------------------------

ccd_transmitance_c1 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 1", "ccd.csv"),
        dtype=np.float64,
        skiprows=1,
    )["(%)"]
    / 100
)


ccd_transmitance_c2 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 2", "ccd.csv"),
        dtype=np.float64,
        skiprows=1,
    )["(%)"]
    / 100
)


ccd_transmitance_c3 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 3", "ccd.csv"),
        dtype=np.float64,
        skiprows=1,
    )["(%)"]
    / 100
)


ccd_transmitance_c4 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 4", "ccd.csv"),
        dtype=np.float64,
        skiprows=1,
    )["(%)"]
    / 100
)

# ------------------------------------ Dichroics ------------------------------------------------------

dichroic_c1_1 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 1", "dichroic_1.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c1_2 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 1", "dichroic_2.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c1_1 = calculate_spline(
    dichroic_c1_1["(%)"], dichroic_c1_1["(nm)"], wavelength_interval
)
dichroic_c1_2 = calculate_spline(
    dichroic_c1_2["(%)"], dichroic_c1_2["(nm)"], wavelength_interval
)
dichroic_c1 = np.multiply(dichroic_c1_1, dichroic_c1_2)

dichroic_c2_1 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 2", "dichroic_1.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c2_2 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 2", "dichroic_2.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c2_1 = calculate_spline(
    dichroic_c2_1["(%)"], dichroic_c2_1["(nm)"], wavelength_interval
)
dichroic_c2_2 = calculate_spline(
    dichroic_c2_2["(%)"], dichroic_c2_2["(nm)"], wavelength_interval
)
dichroic_c2 = np.multiply(dichroic_c2_1, dichroic_c2_2)


dichroic_c3_1 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 3", "dichroic_1.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c3_2 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 3", "dichroic_2.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c3_1 = calculate_spline(
    dichroic_c3_1["(%)"], dichroic_c3_1["(nm)"], wavelength_interval
)
dichroic_c3_2 = calculate_spline(
    dichroic_c3_2["(%)"], dichroic_c3_2["(nm)"], wavelength_interval
)
dichroic_c3 = np.multiply(dichroic_c3_1, dichroic_c3_2)


dichroic_c4_1 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 4", "dichroic_1.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c4_2 = (
    pd.read_csv(
        os.path.join("SPARC4_Spectral_Response", "Channel 4", "dichroic_2.csv"),
        dtype=np.float64,
        skiprows=1,
        decimal=",",
    )
    / 100
)

dichroic_c4_1 = calculate_spline(
    dichroic_c4_1["(%)"], dichroic_c4_1["(nm)"], wavelength_interval
)
dichroic_c4_2 = calculate_spline(
    dichroic_c4_2["(%)"], dichroic_c4_2["(nm)"], wavelength_interval
)
dichroic_c4 = np.multiply(dichroic_c4_1, dichroic_c4_2)


# ------------------------------------ Analysers  ------------------------------------------------------

analyser_ordinary_ray = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "analyser_ordinary.csv"), "rb"),
    delimiter=",",
)

analyser_extra_ordinary_ray = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "analyser_extra_ordinary.csv"), "rb"),
    delimiter=",",
)
