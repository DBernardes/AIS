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
ccd_temp = -70
sky_condition = "photometric"
moon_condition = "new"
wavelength_interval = range(l_init, l_final + l_step, l_step)
wavelength_interval_len = len(wavelength_interval)
sc = Spectrum_Calculation(wavelength_interval=wavelength_interval)
star_specific_photons_per_second = [
    sc.calculate_star_specific_photons_per_second(magnitude)
]
sky_specific_photons_per_second = [
    sc.calculate_sky_specific_photons_per_second(moon_condition)
]
ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "hss": 1,
    "binn": 1,
    "t_exp": 1,
    "ccd_temp": -70,
    "image_size": 1024,
}


sparc4_operation_mode = {"acquisition_mode": "photometric"}


def multiply_matrices(matrix, specific_photons_per_second):
    temp = []
    for array in specific_photons_per_second:
        for idx, value in enumerate(array[0]):
            array[:, idx] = np.dot(matrix, array[:, idx])
        temp.append(array)
    return temp


def calculate_spline(transmitance, component_wavelength_interv, wavelength_interval):
    spl = splrep(component_wavelength_interv, transmitance)
    transmitance = splev(wavelength_interval, spl)
    return transmitance


def convert_magnitude(wavelength, magnitude):
    _H = 6.62607004e-34  # m2 kg / s
    _C = 3e8  # m/s
    _K = 1.38064852e-23  # m2 kg s-2 K-1
    B = 0.2e-6  # m
    S_0 = 4e-2  # W/m2/m
    tel_area = 0.804  # m2
    specific_photons_per_second = (
        S_0 * 10 ** (-magnitude / 2.5) * wavelength * 1e-9 * B * tel_area / (_H * _C)
    )
    return specific_photons_per_second


# ----------------------- importing the moon magnitude ----------------
spreadsheet_path = os.path.join("AIS", "Spectrum_Calculation", "moon_magnitude.csv")
spreadsheet = pd.read_csv(spreadsheet_path, dtype=np.float64)
moon_wavelength_interval = spreadsheet["wavelength"]
moon_magnitude = spreadsheet[moon_condition]
# spl = splrep(moon_wavelength_interval, moon_magnitudes)
# moon_magnitude = splev(moon_wavelength_interval, spl)

# ----------------------- importing the atmosphere spectral response ----------------
spreadsheet_path = os.path.join(
    "AIS", "Atmosphere_Spectral_Response", "atmosphere_spectral_response.csv"
)
spreadsheet = pd.read_csv(spreadsheet_path, dtype=np.float64)
atm_wavelength_interval = spreadsheet["wavelength"]
photometric_extinction_coef = spreadsheet["photometric"] / 100
regular_extinction_coef = spreadsheet["regular"] / 100
good_extinction_coef = spreadsheet["good"] / 100
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
analyser_transmitance = ss["(%)"] / 100
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


ss = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "collimator.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)

collimator_transmitance = (
    calculate_spline(
        ss["(%)"],
        ss["(nm)"],
        wavelength_interval,
    )
    / 100
)


# ------------------------------------ Dichroics ------------------------------------------------------


dichroic_c0 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 0", "dichroic.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


spl = splrep(dichroic_c0["(nm)"], dichroic_c0["(%)"])
dichroic_c0 = splev(wavelength_interval, spl) / 100


dichroic_c1 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 1", "dichroic.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


spl = splrep(dichroic_c1["(nm)"], dichroic_c1["(%)"])
dichroic_c1 = splev(wavelength_interval, spl) / 100


dichroic_c2 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 2", "dichroic.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


spl = splrep(dichroic_c2["(nm)"], dichroic_c2["(%)"])
dichroic_c2 = splev(wavelength_interval, spl) / 100


dichroic_c3 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 3", "dichroic.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


spl = splrep(dichroic_c3["(nm)"], dichroic_c3["(%)"])
dichroic_c3 = splev(wavelength_interval, spl) / 100


dichroic_c4 = pd.read_csv(
    os.path.join("AIS", "SPARC4_Spectral_Response", "Channel 4", "dichroic.csv"),
    dtype=np.float64,
    skiprows=1,
    decimal=".",
)


spl = splrep(dichroic_c4["(nm)"], dichroic_c4["(%)"])
dichroic_c4 = splev(wavelength_interval, spl) / 100

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
