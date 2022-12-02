# Test the sparc4 spectral response
from AIS.Spectral_Response import Channel, Atmosphere, Telescope
import matplotlib.pyplot as plt
import numpy as np
import os
from sys import exit

sr_total = np.ones(200)
wv = np.linspace(400, 1100, 200)
_channel_id = 1
channel = Channel(_channel_id)
sr = channel.get_spectral_response(wv, 'collimator.csv')
sr_total *= sr
plt.plot(wv, sr, label='collimator')
sr = channel.get_spectral_response(wv, f'Channel {_channel_id}/dichroic.csv')
plt.plot(wv, sr, label='dichroic')
sr_total *= sr
sr = channel.get_spectral_response(wv, f'Channel {_channel_id}/camera.csv')
sr_total *= sr
plt.plot(wv, sr, label='camera')
sr = channel.get_spectral_response(wv, f'Channel {_channel_id}/ccd.csv')
sr_total *= sr
plt.plot(wv, sr, label='ccd')

atm = Atmosphere()
sr = atm.get_spectral_response(wv, 1, 'photometric')
sr_total *= sr
plt.plot(wv, sr, label='atmosphere')

tel = Telescope()
sr = tel.get_spectral_response(wv)
sr_total *= sr
plt.plot(wv, sr, label='telescope')


plt.plot(wv, sr_total, 'k', label='total')
plt.legend()
#plt.ylim(0, 1.1)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Transmitance (%)')
plt.title(f'Channel {_channel_id}')
plt.savefig(os.path.join('notebook_figures', f'Channel {_channel_id}.png'))
plt.show()
