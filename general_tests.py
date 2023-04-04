from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

ccd_operation_mode = {
    'em_mode': 'Conv',
    'em_gain': 1,
    'preamp': 1,
    'readout': 1,
    'binn': 1,
    't_exp': 1,
    'image_size': 100
}

ais = Artificial_Image_Simulator(
    ccd_operation_mode, channel_id=1, ccd_temperature=-70)

ais.create_source_sed(calculation_method='blackbody',
                      magnitude=13, wavelength_interval=(400, 1100, 100), temperature=5700)
ais.create_sky_sed(moon_phase='new')
ais.apply_atmosphere_spectral_response(1.5)
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response(
    'photometry')
path = r'C:\Users\observer\Desktop\Testes AIS\data\images'
ais.create_artificial_image(path, (50, 50))
ais.create_bias_image(path)
# ais._integrate_sed()
# print(ais.star_photons_per_second)
