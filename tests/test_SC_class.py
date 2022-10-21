# # -*- coding: utf-8 -*-
# """Spectrum calculation tests.

# This script tests the operation of the spectrum calculation class.
# """

# import numpy as np
# import pytest
# from AIS.Spectrum_Calculation import Spectrum_Calculation

# from .AIS_spectral_response_curves import (
#     calculate_spline,
#     convert_magnitude,
#     magnitude,
#     moon_condition,
#     moon_magnitude,
#     moon_wavelength_interval,
#     wavelength_interval,
#     wavelength_interval_len,
# )

# temperature = 5700


# @pytest.fixture
# def sc():
#     return Spectrum_Calculation(wavelength_interval)


# def test_init(sc):
#     assert sc.wavelength_interval == wavelength_interval


# def test_calculate_specific_photons_per_second(sc):
#     star_specific_photons_per_second = sc.calculate_star_specific_photons_per_second(
#         magnitude
#     )
#     assert np.allclose(
#         star_specific_photons_per_second, star_specific_photons_per_second
#     )


# def test_read_spreadsheet(sc):
#     _, new_moon_magnitude = sc._read_spreadsheet(moon_condition)
#     assert np.allclose(new_moon_magnitude, moon_magnitude)


# def test_calculate_spline(sc):
#     new_moon_magnitudes = sc._calculate_spline(moon_wavelength_interval, moon_magnitude)
#     calculated_moon_magnitudes = calculate_spline(
#         moon_magnitude, moon_wavelength_interval, wavelength_interval
#     )
#     assert np.allclose(new_moon_magnitudes, calculated_moon_magnitudes)


# def test_calculate_sky_specific_photons_per_second(sc):
#     new_moon_magnitudes = calculate_spline(
#         moon_magnitude, moon_wavelength_interval, wavelength_interval
#     )
#     temp = []
#     for magnitude, wavelength in zip(new_moon_magnitudes, wavelength_interval):
#         specific_photons_per_second = convert_magnitude(wavelength, magnitude)
#         temp.append(specific_photons_per_second)

#     sky_specific_photons_per_second = np.zeros((4, wavelength_interval_len))
#     sky_specific_photons_per_second[0, :] = temp

#     new_sky_specific_photons_per_second = sc.calculate_sky_specific_photons_per_second(
#         moon_condition
#     )
#     assert np.allclose(
#         sky_specific_photons_per_second, new_sky_specific_photons_per_second
#     )


# def test_convert_magnitude(sc):
#     _H = 6.62607004e-34  # m2 kg / s
#     _C = 3e8  # m/s
#     _K = 1.38064852e-23  # m2 kg s-2 K-1
#     B = 0.2e-6  # m
#     S_0 = 4e-2  # W/m2/m
#     tel_area = 0.804  # m2
#     Lambda = 350  # nm
#     new_specific_photons_per_second = (
#         S_0 * 10 ** (-magnitude / 2.5) * Lambda * 1e-9 * B * tel_area / (_H * _C)
#     )
#     specific_photons_per_second = sc._convert_magnitude(Lambda, magnitude)

#     assert specific_photons_per_second == new_specific_photons_per_second
