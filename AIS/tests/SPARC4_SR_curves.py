# This file has all the spectral response curves used in the test_SR_SR.py file

import os
from sys import exit

import numpy as np
import pandas as pd
from scipy.interpolate import splev, splrep
from Spectrum_Calculation import Spectrum_Calculation

init, final, step = 400, 1150, 50
magnitude = 22
sc = Spectrum_Calculation(5700, init, final, step)
specific_flux = sc.calculate_specific_flux(magnitude)
wavelength_interval = range(init, final, step)
wavelength_interval_len = len(wavelength_interval)


def calculate_spline(transmitance, component_wavelength_interv, wavelength_interval):
    spl = splrep(component_wavelength_interv, transmitance)
    transmitance = splev(wavelength_interval, spl)
    return transmitance


# ------------------------------------ Calibration Wheel  ------------------------------------------------------

calibration_wheel = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "calibration_wheel.csv"), "rb"),
    delimiter=",",
)

# ------------------------------------ Retarder  ------------------------------------------------------

retarder = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "retarder.csv"), "rb"),
    delimiter=",",
)


# ------------------------------------ Analysers  ------------------------------------------------------

analyser_ordinary_ray = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "analyser_ordinary.csv"), "rb"),
    delimiter=",",
)

analyser_extra_ordinary_ray = np.loadtxt(
    open(os.path.join("SPARC4_Spectral_Response", "analyser_extra_ordinary.csv"), "rb"),
    delimiter=",",
)


# ------------------------------------ Collimator ------------------------------------------------------


colimator_transmitance = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "collimator.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

colimator_transmitance = (
    calculate_spline(
        colimator_transmitance["(%)"],
        colimator_transmitance["(nm)"],
        wavelength_interval,
    )
    / 100
)


# ------------------------------------ Dichroics ------------------------------------------------------


dichroic_c0_1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 0", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c0_2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 0", "dichroic_2.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c0_1 = calculate_spline(
    dichroic_c0_1["(%)"], dichroic_c0_1["(nm)"], wavelength_interval
)
dichroic_c0_2 = calculate_spline(
    dichroic_c0_2["(%)"], dichroic_c0_2["(nm)"], wavelength_interval
)
dichroic_c0 = np.multiply(dichroic_c0_1 / 100, dichroic_c0_2 / 100)


dichroic_c1_1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 1", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c1_2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 1", "dichroic_2.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c1_1 = calculate_spline(
    dichroic_c1_1["(%)"], dichroic_c1_1["(nm)"], wavelength_interval
)
dichroic_c1_2 = calculate_spline(
    dichroic_c1_2["(%)"], dichroic_c1_2["(nm)"], wavelength_interval
)
dichroic_c1 = np.multiply(dichroic_c1_1 / 100, dichroic_c1_2 / 100)

dichroic_c2_1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 2", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c2_2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 2", "dichroic_2.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c2_1 = calculate_spline(
    dichroic_c2_1["(%)"], dichroic_c2_1["(nm)"], wavelength_interval
)
dichroic_c2_2 = calculate_spline(
    dichroic_c2_2["(%)"], dichroic_c2_2["(nm)"], wavelength_interval
)
dichroic_c2 = np.multiply(dichroic_c2_1 / 100, dichroic_c2_2 / 100)


dichroic_c3_1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 3", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c3_2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 3", "dichroic_2.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c3_1 = calculate_spline(
    dichroic_c3_1["(%)"], dichroic_c3_1["(nm)"], wavelength_interval
)
dichroic_c3_2 = calculate_spline(
    dichroic_c3_2["(%)"], dichroic_c3_2["(nm)"], wavelength_interval
)
dichroic_c3 = np.multiply(dichroic_c3_1 / 100, dichroic_c3_2 / 100)


dichroic_c4_1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 4", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c4_2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 4", "dichroic_2.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c4_1 = calculate_spline(
    dichroic_c4_1["(%)"], dichroic_c4_1["(nm)"], wavelength_interval
)
dichroic_c4_2 = calculate_spline(
    dichroic_c4_2["(%)"], dichroic_c4_2["(nm)"], wavelength_interval
)
dichroic_c4 = np.multiply(dichroic_c4_1 / 100, dichroic_c4_2 / 100)


# ------------------------------------ Camera ------------------------------------------------------

camera_c0 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 0", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c0 = (
    calculate_spline(camera_c0["(%)"], camera_c0["(nm)"], wavelength_interval) / 100
)

camera_c1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 1", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c1 = (
    calculate_spline(camera_c1["(%)"], camera_c1["(nm)"], wavelength_interval) / 100
)

camera_c2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 2", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c2 = (
    calculate_spline(camera_c2["(%)"], camera_c2["(nm)"], wavelength_interval) / 100
)


camera_c3 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 4", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c3 = (
    calculate_spline(camera_c3["(%)"], camera_c3["(nm)"], wavelength_interval) / 100
)


camera_c4 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 4", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c4 = (
    calculate_spline(camera_c4["(%)"], camera_c4["(nm)"], wavelength_interval) / 100
)


# ------------------------------------ CCD ------------------------------------------------------

ccd_transmitance_c0 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 0", "ccd.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
ccd_transmitance_c0 = (
    calculate_spline(
        ccd_transmitance_c0["(%)"], ccd_transmitance_c0["(nm)"], wavelength_interval
    )
    / 100
)


ccd_transmitance_c1 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 1", "ccd.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


ccd_transmitance_c1 = (
    calculate_spline(
        ccd_transmitance_c1["(%)"], ccd_transmitance_c1["(nm)"], wavelength_interval
    )
    / 100
)


ccd_transmitance_c2 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 2", "ccd.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

ccd_transmitance_c2 = (
    calculate_spline(
        ccd_transmitance_c2["(%)"], ccd_transmitance_c2["(nm)"], wavelength_interval
    )
    / 100
)


ccd_transmitance_c3 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 3", "ccd.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

ccd_transmitance_c3 = (
    calculate_spline(
        ccd_transmitance_c3["(%)"], ccd_transmitance_c3["(nm)"], wavelength_interval
    )
    / 100
)


ccd_transmitance_c4 = pd.read_csv(
    os.path.join("SPARC4_Spectral_Response", "Channel 4", "ccd.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

ccd_transmitance_c4 = (
    calculate_spline(
        ccd_transmitance_c4["(%)"], ccd_transmitance_c4["(nm)"], wavelength_interval
    )
    / 100
)
