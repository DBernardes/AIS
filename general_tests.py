# -*- coding: utf-8 -*-

"""
Created on Tue Apr 27 10:23:27 2021

This is a script for the development of general tests

@author: denis
"""

# -----------------------------------------------------------------------------------------
# Test the PSF in hte polarimetric mode
# from AIS.Point_Spread_Function import Point_Spread_Function
# from tests.AIS_spectral_response_curves import ccd_operation_mode
# import matplotlib.pyplot as plt
# psf = Point_Spread_Function(ccd_operation_mode, 1)
# image = psf.create_star_image((50, 50), 100, 100)
# plt.imshow(image)
# plt.show()
# -----------------------------------------------------------------------------------------
# from AIS.Header import Header
# Pint the header
# from tests.AIS_spectral_response_curves import ccd_operation_mode
# hdr = Header(ccd_operation_mode)
# header = hdr.create_header()
# print(repr(header))
# -----------------------------------------------------------------------------------------
# Test the blackbody profile

# from AIS. Artificial_Image_Simulator import Artificial_Image_Simulator
# from AIS.Spectral_Energy_Distribution import Source
# import matplotlib.pyplot as plt
# from tests.AIS_spectral_response_curves import ccd_operation_mode

# ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)

# ais.create_source_sed('blackbody', 15, (400, 1100, 100), 5700)
# ais.create_sky_sed('new', ais.wavelength)
# plt.plot(ais.wavelength, ais.source_sed, label='Blackbody')
# ais.apply_atmosphere_spectral_response()
# plt.plot(ais.wavelength, ais.source_sed, label='Blackbody with atmosphere')
# ais.apply_telescope_spectral_response()
# plt.plot(ais.wavelength, ais.source_sed,
#          label='Blackbody with atmosphere and telescope')
# ais.apply_sparc4_spectral_response('photometric')
# plt.plot(ais.wavelength, ais.source_sed,
#          label='Blackbody with atmosphere, telescope and SPARC4')
# plt.legend()
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('SED (W/m2/nm)')
# plt.title('Blackbody SED')
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the sbpy package
# from sbpy.calib import Vega, vega_fluxd

# print(vega_fluxd.get()["Johnson V"].value)
# -----------------------------------------------------------------------------------------
# Test the _read_spectral_library function
# from AIS.Spectral_Energy_Distribution import Source
# from matplotlib import pyplot as plt
# source = Source()
# wv, sed = source.calculate_sed(
#     'spectral_library', 10, spectral_type='G')
# plt.plot(wv, sed)
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the Sky Class
# from AIS.Spectral_Energy_Distribution import Sky
# import matplotlib.pyplot as plt
# from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength
# sky = Sky()
# sed = sky.calculate_sed('new', obj_wavelength)
# plt.plot(obj_wavelength, sed)
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the sparc4 spectral response
# from AIS.Spectral_Response import Channel
# import matplotlib.pyplot as plt
# import numpy as np

# sed = np.ones(100)
# wv = np.linspace(400, 1100, 100)
# for i in [1, 2, 3, 4]:
#     channel = Channel(i)
#     channel.write_sparc4_operation_mode('photometric')
#     ch1_spectral_response = channel.apply_spectral_response(sed, wv)
#     plt.plot(wv, ch1_spectral_response, label=f'Channel {i}')
# plt.legend()
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Transmitance (%)')
# plt.show()
# -----------------------------------------------------------------------------------------
# Test the atmosphere spectral response
# from AIS.Spectral_Response import Atmosphere
# import matplotlib.pyplot as plt
# from tests.AIS_spectral_response_curves import wavelength_interval as obj_wavelength

# atmosphere = Atmosphere()
# spectral_response = atmosphere.get_spectral_response(obj_wavelength, 1)
# plt.plot(obj_wavelength, spectral_response)
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Transmitance (%)')
# plt.show()
# -----------------------------------------------------------------------------------------
# Test telescope spectral response

# telescope = Telescope()
# spectral_response = telescope.get_spectral_response(obj_wavelength)
# plt.plot(obj_wavelength, spectral_response)
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Transmitance (%)')
# plt.show()


# -----------------------------------------------------------------------------------------
# Test the AIS.create_image function
from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
from tests.AIS_spectral_response_curves import ccd_operation_mode
import matplotlib.pyplot as plt

wavelegnth_interval = (350, 1100, 100)
ais = Artificial_Image_Simulator(ccd_operation_mode, 1, -70)
ais.create_source_sed('blackbody', 30, wavelegnth_interval, 5700)
ais.create_sky_sed('new')
ais.apply_atmosphere_spectral_response()
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response('photometric')
ais._integrate_sed()
print(ais.star_photons_per_second)
#ais.create_artificial_image(r'E:\images\test', (50, 50))
