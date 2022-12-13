from AIS.Spectral_Response import Channel
import numpy as np
import matplotlib.pyplot as plt
from sys import exit

n = 100
wv = np.linspace(300, 1000, n)
sed = np.ones(n)
ch = Channel(1)
ch.write_sparc4_operation_mode(acquisition_mode='polarimetry', calibration_wheel='polarizer',
                               retarder_waveplate='half', retarder_waveplate_angle=30)

sed = ch.apply_spectral_response(sed, wv)
plt.plot(wv, sed[0], 'b')
plt.plot(wv, sed[1], 'k')
plt.show()
