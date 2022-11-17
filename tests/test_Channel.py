# # -*- coding: utf-8 -*-
# Tests of the Channel Class
# Oct 17th 2022
# @author: denis
#

import os
import numpy as np
import pandas as pd
import pytest
from AIS.Spectral_Response import Channel
from scipy.interpolate import splev, splrep

from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength

CHANNEL = 1
BASE_PATH = os.path.join("AIS", "Spectral_Response", "channel")
_POL_OPTICAL_COMPONENTS = {
    "polarizer": "polarizer.csv",
    "depolarizer": "depolarizer.csv",
    "retarder": "retarder.csv",
    "analyser": "analyser.csv",
}

_PHOT_OPTICAL_COMPONENTS = {
    "collimator": "collimator.csv",
    "dichroic": "dichroic.csv",
    "camera": "camera.csv",
    "ccd": "ccd.csv",
}


@pytest.fixture
def channel():
    return Channel(CHANNEL, 'photometric')


def test_csv_file_name_pol(channel):
    assert channel._POL_OPTICAL_COMPONENTS == _POL_OPTICAL_COMPONENTS


def test_csv_file_name_phot(channel):
    assert channel._PHOT_OPTICAL_COMPONENTS == _PHOT_OPTICAL_COMPONENTS


def test_base_path(channel):
    assert channel._BASE_PATH == BASE_PATH


# --------------------------------------------------------------------------------------------

polarizer_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["polarizer"])
ss = pd.read_csv(polarizer_path)
wv_polarizer = ss["Wavelength (nm)"]
sr_polarizer = ss["Transmitance (%)"] / 100


def test_read_csv_file_polarizer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(polarizer_path)
    assert np.allclose(new_wavelength, wv_polarizer)
    assert np.allclose(new_transmitance, sr_polarizer)


depolarizer_path = os.path.join(
    BASE_PATH, _POL_OPTICAL_COMPONENTS["depolarizer"])
ss = pd.read_csv(depolarizer_path)
wv_depolarizer = ss["Wavelength (nm)"]
sr_depolarizer = ss["Transmitance (%)"] / 100


def test_read_csv_file_depolarizer(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(depolarizer_path)
    assert np.allclose(new_wavelength, wv_depolarizer)
    assert np.allclose(new_transmitance, sr_depolarizer)


retarder_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["retarder"])
ss = pd.read_csv(retarder_path)
wv_retarder = ss["Wavelength (nm)"]
sr_retarder = ss["Transmitance (%)"] / 100


def test_read_csv_file_retarder(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(retarder_path)
    assert np.allclose(new_wavelength, wv_retarder)
    assert np.allclose(new_transmitance, sr_retarder)


analyser_path = os.path.join(BASE_PATH, _POL_OPTICAL_COMPONENTS["analyser"])
ss = pd.read_csv(analyser_path)
wv_analyser = ss["Wavelength (nm)"]
sr_analyser = ss["Transmitance (%)"] / 100


def test_read_csv_file_analyser(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(analyser_path)
    assert np.allclose(new_wavelength, wv_analyser)
    assert np.allclose(new_transmitance, sr_analyser)


collimator_path = os.path.join(
    BASE_PATH, _PHOT_OPTICAL_COMPONENTS["collimator"])
ss = pd.read_csv(collimator_path)
wv_collimator = ss["Wavelength (nm)"]
sr_collimator = ss["Transmitance (%)"] / 100


def test_read_csv_file_collimator(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(collimator_path)
    assert np.allclose(new_wavelength, wv_collimator)
    assert np.allclose(new_transmitance, sr_collimator)


dichroic_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["dichroic"]
)
ss = pd.read_csv(dichroic_path)
wv_dichroic = ss["Wavelength (nm)"]
sr_dichroic = ss["Transmitance (%)"] / 100


def test_read_csv_file_dichroic(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(dichroic_path)
    assert np.allclose(new_wavelength, wv_dichroic)
    assert np.allclose(new_transmitance, sr_dichroic)


camera_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["camera"])
ss = pd.read_csv(camera_path)
wv_camera = ss["Wavelength (nm)"]
sr_camera = ss["Transmitance (%)"] / 100


def test_read_csv_file_camera(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(camera_path)
    assert np.allclose(new_wavelength, wv_camera)
    assert np.allclose(new_transmitance, sr_camera)


ccd_path = os.path.join(
    BASE_PATH, f"Channel {CHANNEL}", _PHOT_OPTICAL_COMPONENTS["ccd"])
ss = pd.read_csv(ccd_path)
wv_ccd = ss["Wavelength (nm)"]
sr_ccd = ss["Transmitance (%)"] / 100


def test_read_csv_file_ccd(channel):
    new_wavelength, new_transmitance = channel._read_csv_file(ccd_path)
    assert np.allclose(new_wavelength, wv_ccd)
    assert np.allclose(new_transmitance, sr_ccd)

# -----------------------Test SPARC4 operation mode ---------------------------------------


def test_verify_sparc4_operation_mode_photometric():
    channel = Channel(CHANNEL, 'photometric')
    channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_error():
    channel = Channel(CHANNEL, 'error')
    with pytest.raises(ValueError):
        channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric():
    channel = Channel(CHANNEL, "polarimetric", 'empty', 'half')
    channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric_polarizer():
    channel = Channel(CHANNEL, "polarimetric", 'polarizer', 'half')
    channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric_depolarizer():
    channel = Channel(CHANNEL, "polarimetric", 'depolarizer', 'half')
    channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric_quarter():
    channel = Channel(CHANNEL, "polarimetric", 'empty', 'quarter')
    channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric_wrong_retarder_value():
    channel = Channel(CHANNEL, "polarimetric", 'empty', '')
    with pytest.raises(ValueError):
        channel._verify_sparc4_operation_mode()


def test_verify_sparc4_operation_mode_polarimetric_wrong_calibration_wheel_value():
    channel = Channel(CHANNEL, "polarimetric", '', 'half')
    with pytest.raises(ValueError):
        channel._verify_sparc4_operation_mode()


# --------------------------------------------------------------------------------------------


spl = splrep(wv_collimator, sr_collimator)
new_spectral_response = splev(obj_wavelength, spl)


def test_interpolate_spectral_response(channel):
    class_spectral_response = channel._interpolate_spectral_response(
        wv_collimator, sr_collimator, obj_wavelength
    )
    assert np.allclose(class_spectral_response, new_spectral_response)


def test_get_spectral_response(channel):
    class_spectral_response = channel.get_spectral_response(
        obj_wavelength, _PHOT_OPTICAL_COMPONENTS["collimator"]
    )
    assert np.allclose(class_spectral_response,
                       new_spectral_response, rtol=0.005)


def test_apply_phot_spectral_response(channel):
    n = len(obj_wavelength)
    esd = np.linspace(100, 360, n)
    spl = splrep(wv_collimator, sr_collimator)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(esd, new_spectral_response)
    spl = splrep(wv_dichroic, sr_dichroic)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_camera, sr_camera)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_ccd, sr_ccd)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    class_reduced_esd = channel._apply_photometric_spectral_response(
        esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)


def test_apply_pol_spectral_response_1(channel):
    n = len(obj_wavelength)
    esd = np.linspace(100, 360, n)
    spl = splrep(wv_depolarizer, sr_depolarizer)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(esd, new_spectral_response)
    spl = splrep(wv_retarder, sr_retarder)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_analyser, sr_analyser)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    channel.polarimetric_config['calibration_wheel'] = 'depolarizer'
    class_reduced_esd = channel._apply_polarimetric_spectral_response(
        esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)


def test_apply_pol_spectral_response_2(channel):
    n = len(obj_wavelength)
    esd = np.linspace(100, 360, n)
    spl = splrep(wv_polarizer, sr_polarizer)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(esd, new_spectral_response)
    spl = splrep(wv_retarder, sr_retarder)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_analyser, sr_analyser)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    channel.polarimetric_config['calibration_wheel'] = 'polarizer'
    class_reduced_esd = channel._apply_polarimetric_spectral_response(
        esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)


def test_apply_pol_spectral_response_3(channel):
    n = len(obj_wavelength)
    esd = np.linspace(100, 360, n)
    spl = splrep(wv_retarder, sr_retarder)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(esd, new_spectral_response)
    spl = splrep(wv_analyser, sr_analyser)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    channel.polarimetric_config['calibration_wheel'] = 'empty'
    class_reduced_esd = channel._apply_polarimetric_spectral_response(
        esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)


def test_apply_spectral_response(channel):
    n = len(obj_wavelength)
    esd = np.linspace(100, 360, n)
    spl = splrep(wv_retarder, sr_retarder)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(esd, new_spectral_response)
    spl = splrep(wv_analyser, sr_analyser)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_collimator, sr_collimator)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_dichroic, sr_dichroic)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_camera, sr_camera)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    spl = splrep(wv_ccd, sr_ccd)
    new_spectral_response = splev(obj_wavelength, spl)
    reduced_esd = np.multiply(reduced_esd, new_spectral_response)
    channel.polarimetric_config['calibration_wheel'] = 'empty'
    channel.acquisition_mode = 'polarimetric'
    class_reduced_esd = channel.apply_spectral_response(
        esd, obj_wavelength)
    assert np.allclose(class_reduced_esd, reduced_esd)
