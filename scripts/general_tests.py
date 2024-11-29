import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import h, c
import pandas as pd
from scipy.optimize import curve_fit
from math import atan2, sqrt


ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}
_TELESCOPE_EFFECTIVE_AREA = 0.804  # m2


def convert_from_flux_to_photons_per_sec(csv_file):
    spec = pd.read_csv(csv_file)
    source_wavelength = np.asarray(spec["wavelength"])
    source_sed = np.asarray(spec["flux"]) * 1e9
    for idx, _ in enumerate(source_sed):
        source_sed[idx] *= (
            _TELESCOPE_EFFECTIVE_AREA
            * source_wavelength[idx]
            * 1e-9
            / (h * c)  # photons/m/s
        )

    return source_wavelength, source_sed


def half_model(psi, q, u):
    return q * np.cos(4 * psi) + u * np.sin(4 * psi)


def calc_zi(ord, extra, k):
    zi = (ord - extra * k) / (extra * k + ord)
    return zi


def calc_k(ord, extra):
    return sum(extra[:4]) / sum(ord[:4])


def apply_half_model(ss):
    angles, ord, extra = np.asarray(ss["angle"]), ss["ord"], ss["extra"]

    k = calc_k(ord, extra)
    zi = calc_zi(ord, extra, k)
    (q, u), _ = curve_fit(half_model, np.deg2rad(angles), zi, method="trf")

    p = sqrt(q**2 + u**2)
    theta = np.rad2deg(atan2(u, q) / 2)
    return p, theta, q, u, k


from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

csv_folder = os.path.join("..", "..", "AIS-tests", "polarimetric tests", "stars", "csv")

csv_file = os.path.join(csv_folder, "stars_info.csv")
stars_info = pd.read_csv(csv_file)
csv_file = os.path.join(csv_folder, "stars", "Vela 1-95", "gaia spectrum.csv")
wv, sed = convert_from_flux_to_photons_per_sec(csv_file)

vela = stars_info[stars_info["star"] == "Vela 1-95"].reset_index()
BVRI = ["B", "V", "R", "I"]
pBVRI = {k: vela.loc[0, f"pol {k}"] for k in BVRI}

ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
dict = {"angle": [], "ord": [], "extra": []}
for i in range(16):
    # ais.write_source_sed(wv.copy(), sed.copy())
    # ais.apply_Serkowski_curve(pBVRI, 0)
    ais.create_source_sed_blackbody(12, (350, 1100, 100), 5700)
    ais.apply_linear_polarization(1, 20)
    ais.create_sky_sed("new")
    ais.apply_sparc4_spectral_response("polarimetry", "", "ideal-half", i * 22.5)
    ais._integrate_sed()

    dict["angle"].append(i * 22.5)
    dict["ord"].append(ais.star_photons_per_second[0])
    dict["extra"].append(ais.star_photons_per_second[1])


ss = pd.DataFrame.from_dict(dict)

a = apply_half_model(ss)
print(a)
