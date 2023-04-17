# Test polarization with a non-ideal polarizer and waveplate
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
for i in range(16):
    ais.create_source_sed(calculation_method='blackbody',
                          magnitude=15, wavelength_interval=(400, 1100, 100), temperature=5700)
    ais.create_sky_sed(moon_phase='new')
    ais.apply_atmosphere_spectral_response()
    ais.apply_telescope_spectral_response()
    ais.apply_sparc4_spectral_response(
        'polarimetry', 'polarizer', 'quarter', retarder_waveplate_angle=i*22.5)
    ais._integrate_sed()
    print(i*22.5, *ais.star_photons_per_second, sep=',')
