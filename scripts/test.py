import os
import pandas as pd
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

BVRI = ["B", "V", "R", "I"]
ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}

csv_path = os.path.join(
    "..", "..", "AIS-tests", "polarimetric tests", "data", "stars", "csv"
)
# csv with the stars polarization in the BVRI
csv_file = os.path.join(csv_path, "stars_info.csv")
star_info = pd.read_csv(csv_file)

# Polarimetric measurements with SPARC4
csv_file = os.path.join(csv_path, "observations.csv")
observations = pd.read_csv(csv_file)
observations["calwheel"] = observations["calwheel"].replace("none", "")
row = observations.iloc[28]


star = row["star"]
waveplate = row["waveplate"]
channel = row["channel"]
calibration_wheel = row["calwheel"]
stars_info_row = star_info[star_info["star"] == star]


pol_BVRI = {k: stars_info_row[f"pol {k}"].values[0] for k in BVRI}
dict = {"angle": [], "ord": [], "extra": []}


ais = Artificial_Image_Simulator(
    ccd_operation_mode, channel_id=channel, ccd_temperature=-70
)
ais.create_source_sed_spectral_library(
    magnitude=12,
    wavelength_interval=(350, 1100, 100),
    spectral_type=stars_info_row["chosen sp"].values[0],
)
ais.apply_Serkowski_curve(pol_BVRI)
ais.create_sky_sed(moon_phase="new")
ais.apply_atmosphere_spectral_response(row["airmass"], row["sky condition"])
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response(
    "polarimetry",
    calibration_wheel=calibration_wheel,
    retarder_waveplate=waveplate,
    retarder_waveplate_angle=1 * 22.5,
)
