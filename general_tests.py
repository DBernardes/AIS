from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator

# Create an instance of the Artificial_Image_Simulator class

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
ais.create_source_sed('blackbody', 15, (400, 1100, 100), 5700)
ais.create_sky_sed('new')
ais.apply_atmosphere_spectral_response(
    air_mass=1.2, sky_condition='photometric')
ais.apply_telescope_spectral_response()
ais.apply_sparc4_spectral_response('photometry')
ais.create_background_image(r'E:\images\test')
