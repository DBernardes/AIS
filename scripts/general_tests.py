# APPLY SERKOWSKI POLARIZATION CURVE

import os
from math import atan2, pi, sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Spectral_Response._utils import calculate_polarizer_matrix

path = os.path.join("..", "AIS-tests", "polarimetric tests", "data", "stars")
ss = pd.read_csv(os.path.join(path, "stars_info.csv"))
pol_values = {"B": 8.602, "V": 9.548, "R": 9.671, "I": 9.009}


def half_model(psi, q, u):
    return q * np.cos(4 * psi) + u * np.sin(4 * psi)


def calc_zi(ord, extra, k):
    zi = (ord - extra * k) / (ord + extra * k)
    return zi


def apply_half_model(angles, ord, extra):

    zi = calc_zi(ord, extra, 1)

    (q, u), pcov = curve_fit(half_model, np.deg2rad(angles), zi, method="trf")
    p = sqrt(q**2 + u**2)
    theta = np.rad2deg(atan2(u, q)) / 2
    if theta < 0:
        theta += 180
    return q, u, p, theta


ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 60,
    "image_size": 100,
}

dict = {"angle": [], "ord": [], "extra": []}
ais = Artificial_Image_Simulator(ccd_operation_mode, channel_id=1, ccd_temperature=-70)
ais.create_source_sed(
    calculation_method="spectral_library",
    magnitude=12.17,
    wavelength_interval=(463.9, 464, 2),
    spectral_type="b0v",
)
ais.apply_Serkowski_curve(pol_values)
ais.create_sky_sed(moon_phase="new")
ais.apply_atmosphere_spectral_response(1.101, "photometric")
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response(
    "polarimetry",
    retarder_waveplate="ideal-half",
    retarder_waveplate_angle=0,
)
ais._integrate_sed()
ord = ais.star_photons_per_second[0]
extra = ais.star_photons_per_second[1]
print(ord, extra)
print(calc_zi(ord, extra, 1))

# * expecting 9.08
# * Como a medida é uma média, deve ficar abaixo?
# * 90862.69469041278 76621.16026482507
# * 0.0850322822423328

# angles, ord, extra = (
#     np.asarray(dict["angle"]),
#     np.asarray(dict["ord"]),
#     np.asarray(dict["extra"]),
# )
# q, u, p, theta = apply_half_model(angles, ord, extra)


# plt.title(f"q={q:.3f}; u={u:.3f}; p={p:.2%}")
# zi = calc_zi(ord, extra, 1)
# plt.plot(angles, zi, "bo", label="dados AIS")
# angles = np.linspace(0, 2 * pi, 100)
# model = half_model(angles, q, u)
# plt.plot(np.rad2deg(angles), model, "b-", label="modelo")
# plt.xlabel("Ângulo de lâmina retardadora (graus)")
# plt.ylabel("Amplitude da modulação")
# plt.legend()
# plt.show()
