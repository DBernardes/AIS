import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from AIS.Spectral_Energy_Distribution import Source
from AIS.Spectral_Response import Channel

ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}


# ais = Artificial_Image_Simulator(ccd_operation_mode, channel_id=1, ccd_temperature=-70)


# ais.create_source_sed(
#     calculation_method="spectral_library",
#     magnitude=15,
#     wavelength_interval=(400, 1100, 100),
#     spectral_type="m0iii",
# )
# ais.apply_linear_polarization(100)
# ais.create_sky_sed(moon_phase="new")
# ais.apply_sparc4_spectral_response(
#     "polarimetry",
#     retarder_waveplate="quarter",
#     retarder_waveplate_angle=22.5,
# )
# ais._integrate_sed()
# ord, extra = ais.star_photons_per_second[0], ais.star_photons_per_second[1]
# print(22.5, ord, extra)


sed = Source()
wv, _ = sed.calculate_sed(
    calculation_method="spectral_library",
    magnitude=15,
    wavelength_interval=(400, 1100, 100),
    spectral_type="m0iii",
)
sed = sed.apply_linear_polarization(100)

ch = Channel(1)
ch.write_sparc4_operation_mode(
    "polarimetry",
    retarder_waveplate="ideal-quarter",
    retarder_waveplate_angle=20,
)

ch.sed = sed
ch.obj_wavelength = wv
ch._apply_polarimetric_spectral_response()
ch._apply_photometric_spectral_response()
# sed = ch.sed
# ord = np.trapz(sed[0], wv * 1e-9)
# extra = np.trapz(sed[1], wv * 1e-9)
# print((ord - extra) / (ord + extra))


# * Investigando o caso de 22.5°
# * Acontece para um caso específico do espectro (canais 1 e 2)
# * Fixando diferença de fase em 90°, tudo certo
# * Fazendo a diferença de fase 89° aparece o problema
# * Não seria o mesmo caso da L4 ideal para 20° ?
