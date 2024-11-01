# PLOT SPECTRAL LIBRARY WITH ATM + TELESCOPE + SPARC4
import matplotlib.pyplot as plt

from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

ccd_operation_mode = {
    "em_mode": "Conv",
    "em_gain": 1,
    "preamp": 1,
    "readout": 1,
    "binn": 1,
    "t_exp": 1,
    "image_size": 100,
}

ais = Artificial_Image_Simulator(ccd_operation_mode, 2, -70)

ais.create_source_sed_spectral_library(15, (350, 1100, 1000), spectral_type="g5v")
ais.create_sky_sed("new")
plt.plot(ais.wavelength, ais.source_sed[0], label="G5 v")
ais.apply_atmosphere_spectral_response()
plt.plot(ais.wavelength, ais.source_sed[0], label="G5 v + atm")
ais.apply_telescope_spectral_response()
plt.plot(ais.wavelength, ais.source_sed[0], label="G5 v + atm + tel")
ais.apply_sparc4_spectral_response("photometry")
plt.plot(ais.wavelength, ais.source_sed, label="G5 v + atm + tel + inst")
ais._integrate_sed()
print(f"{ais.star_photons_per_second:.2f}")
plt.legend()
plt.xlabel("Wavelength (nm)")
plt.ylabel(r"SED (photons/m$^2$/nm)")
plt.xlim(350, 1100)
plt.ylim(0, 1e11)
plt.savefig("notebook/figures/spectral_library.png", dpi=300)
plt.show()
