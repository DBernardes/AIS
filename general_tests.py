from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

ccd_operation_mode = {
    'em_mode': 'Conv', 'em_gain': 1, 'preamp': 1, 'readout': 1, 'binn': 1, 't_exp': 1, 'image_size': 100
}

ais = Artificial_Image_Simulator(ccd_operation_mode, channel_id=4, ccd_temperature=-70)
ais.create_source_sed(calculation_method='spectral_library',
                        magnitude=14, wavelength_interval=(400, 1100, 100), spectral_type='G0v')
ais.create_sky_sed(moon_phase='new')
ais.apply_atmosphere_spectral_response()
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response(acquisition_mode='photometry')
t = ais.calculate_exposure_time()
print(t)