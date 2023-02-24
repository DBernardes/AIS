from AIS.Artificial_Image_Simulator import Artificial_Image_Simulator
import matplotlib.pyplot as plt

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
                      magnitude=12, wavelength_interval=(400, 1100, 100), temperature=5700)

plt.plot(ais.wavelength, ais.source_sed)
plt.show()
