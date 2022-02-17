"""This file has all the spectral response curves used in the test_SR_SR.py file"""

import os

import numpy as np
import pandas as pd
from AIS.Spectrum_Calculation import Spectrum_Calculation
from numpy import cos, pi, sin
from scipy.interpolate import splev, splrep

THETA_POL = np.deg2rad(0)
l_init, l_final, l_step = 400, 1100, 50
magnitude = 22
air_mass = 1
sky_condition = "photometric"
wavelength_interval = range(l_init, l_final + l_step, l_step)
wavelength_interval_len = len(wavelength_interval)
sc = Spectrum_Calculation(
    wavelength_interval=wavelength_interval, star_temperature=5700
)
star_specific_photons_per_second = [sc.calculate_specific_photons_per_second(magnitude)]
sky_specific_photons_per_second = [
    sc.calculate_specific_photons_per_second(magnitude + 3)
]


def calculate_spline(transmitance, component_wavelength_interv, wavelength_interval):
    spl = splrep(component_wavelength_interv, transmitance)
    transmitance = splev(wavelength_interval, spl)
    return transmitance


# ----------------------- importing the atmosphere spectral response ----------------
spreadsheet_path = os.path.join(
    "AIS", "Atmosphere_Spectral_Response", "atmosphere_spectral_response.csv"
)
spreadsheet = pd.read_csv(spreadsheet_path)
atm_wavelength_interval = [float(value) for value in spreadsheet["Wavelength"][1:]]
photometric_extinction_coef = [
    float(value) / 100 for value in spreadsheet["photometric"][1:]
]
regular_extinction_coef = [float(value) / 100 for value in spreadsheet["regular"][1:]]
good_extinction_coef = [float(value) / 100 for value in spreadsheet["good"][1:]]
transmitance = [10 ** (-0.4 * k * air_mass) for k in photometric_extinction_coef]
spl = splrep(atm_wavelength_interval, transmitance)
atm_transmitance = splev(wavelength_interval, spl)


# ----------------------- importing the telescope spectral response ----------------
ss = pd.read_csv(
    os.path.join(
        "AIS", "Telescope_Spectral_Response", "telescope_spectral_response.csv"
    ),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

tel_wavelength_interval = ss["(nm)"]
tel_reflectance_spreadsheet = ss["(%)"] / 100
spl = splrep(tel_wavelength_interval, tel_reflectance_spreadsheet)
tel_reflectance = splev(wavelength_interval, spl)

# ------------------------------------ Calibration Wheel  ------------------------------------------------------

ss = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "polarizer.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
polarizer_wavelength_interval = ss["(nm)"]
polarizer_transmitance = ss["(%)"] / 100
spl = splrep(polarizer_wavelength_interval, polarizer_transmitance)
polarizer_transmitance = splev(wavelength_interval, spl)


ss = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "depolarizer.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
depolarizer_wavelength_interval = ss["(nm)"]
depolarizer_transmitance = ss["(%)"] / 100
spl = splrep(depolarizer_wavelength_interval, depolarizer_transmitance)
depolarizer_transmitance = splev(wavelength_interval, spl)


POLARIZER_MATRIX = 0.5 * np.asarray(
    [
        [1, cos(2 * THETA_POL), sin(2 * THETA_POL), 0],
        [
            cos(2 * THETA_POL),
            cos(2 * THETA_POL) ** 2,
            cos(2 * THETA_POL) * sin(2 * THETA_POL),
            0,
        ],
        [
            sin(2 * THETA_POL),
            cos(2 * THETA_POL) * sin(2 * THETA_POL),
            sin(2 * THETA_POL) ** 2,
            0,
        ],
        [0, 0, 0, 0],
    ]
)


# ------------------------------------ Retarder  ------------------------------------------------------


def calculate_retarder_matrix(phase_difference, THETA_POL):
    phase_difference = np.deg2rad(phase_difference)
    retarder_matrix = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * THETA_POL) ** 2
                + sin(2 * THETA_POL) ** 2 * cos(phase_difference),
                cos(2 * THETA_POL) * sin(2 * THETA_POL) * (1 - cos(phase_difference)),
                -sin(2 * THETA_POL) * sin(phase_difference),
            ],
            [
                0,
                cos(2 * THETA_POL) + sin(2 * THETA_POL) * (1 - cos(phase_difference)),
                sin(2 * THETA_POL) ** 2
                + cos(2 * THETA_POL) ** 2 * cos(phase_difference),
                cos(2 * THETA_POL) * sin(phase_difference),
            ],
            [
                0,
                sin(2 * THETA_POL) * sin(phase_difference),
                -cos(2 * THETA_POL) * sin(phase_difference),
                cos(phase_difference),
            ],
        ],
    )

    return retarder_matrix


ss = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "retarder_transmitance.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
retarder_wavelength_interval = ss["(nm)"]
retarder_transmitance = ss["(%)"] / 100
spl = splrep(retarder_wavelength_interval, retarder_transmitance)
retarder_transmitance = splev(wavelength_interval, spl)

file = os.path.join(
    "AIS", "SPARC4_Spectral_Response", "retarder_phase_diff_" + "quarter" + ".csv"
)
ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
retarder_wavelength = ss["(nm)"]
retardance = ss["(waves)"]
spl = splrep(retarder_wavelength, retardance)
retardance = splev(wavelength_interval, spl)
retardance_quarter = [val * 360 for val in retardance]

file = os.path.join(
    "AIS", "SPARC4_Spectral_Response", "retarder_phase_diff_" + "half" + ".csv"
)
ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
retarder_wavelength = ss["(nm)"]
retardance = ss["(waves)"]
spl = splrep(retarder_wavelength, retardance)
retardance = splev(wavelength_interval, spl)
retardance_half = [val * 360 for val in retardance]

# ------------------------------------ Analysers  ------------------------------------------------------

file = os.path.join("AIS", "SPARC4_Spectral_Response", "analyser.csv")
ss = pd.read_csv(file, dtype=np.float64, skiprows=1, decimal=".")
analyser_wavelength = ss["(nm)"]
analyser_transmitance = ss["(%)"]
spl = splrep(analyser_wavelength, analyser_transmitance)
analyser_transmitance = splev(wavelength_interval, spl)

THETA_POL = np.deg2rad(THETA_POL + pi / 2)
POLARIZER_90_MATRIX = 0.5 * np.asarray(
    [
        [1, cos(2 * THETA_POL), sin(2 * THETA_POL), 0],
        [
            cos(2 * THETA_POL),
            cos(2 * THETA_POL) ** 2,
            cos(2 * THETA_POL) * sin(2 * THETA_POL),
            0,
        ],
        [
            sin(2 * THETA_POL),
            cos(2 * THETA_POL) * sin(2 * THETA_POL),
            sin(2 * THETA_POL) ** 2,
            0,
        ],
        [0, 0, 0, 0],
    ]
)


# ------------------------------------ Collimator ------------------------------------------------------


colimator_transmitance = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "collimator.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 0", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c0_2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 0", "dichroic_2.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 1", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c1_2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 1", "dichroic_2.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 2", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c2_2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 2", "dichroic_2.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 3", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c3_2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 3", "dichroic_2.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "dichroic_1.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

dichroic_c4_2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "dichroic_2.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 0", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c0 = (
    calculate_spline(camera_c0["(%)"], camera_c0["(nm)"], wavelength_interval) / 100
)

camera_c1 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 1", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c1 = (
    calculate_spline(camera_c1["(%)"], camera_c1["(nm)"], wavelength_interval) / 100
)

camera_c2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 2", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c2 = (
    calculate_spline(camera_c2["(%)"], camera_c2["(nm)"], wavelength_interval) / 100
)


camera_c3 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c3 = (
    calculate_spline(camera_c3["(%)"], camera_c3["(nm)"], wavelength_interval) / 100
)


camera_c4 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "camera.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)
camera_c4 = (
    calculate_spline(camera_c4["(%)"], camera_c4["(nm)"], wavelength_interval) / 100
)


# ------------------------------------ CCD ------------------------------------------------------

ccd_transmitance_c0 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 0", "ccd.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 1", "ccd.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 2", "ccd.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 3", "ccd.csv"),
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
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "ccd.csv"),
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
