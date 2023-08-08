"""
Noise Class
===========

This class calculates the read noise and the dark noise of the SPARC4 EMCCDs as a function of their operation mode. 
The calculations are done based on a series of characterization of the noise of the these cameras. 
For the conventional mode, the values of the read noise in the Tabelas_Valores_Ruido_Leitura spreadsheet are used. 
For the EM mode, an interpolation of the data presented by the respective spreadshhet is done, as a function of the EM gain.
"""

# Denis Varise Bernardes.
# 08/10/2019.

import os

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

__all__ = ["Noise"]


class Noise:
    """Noise Class.

    This class calculates the read noise and the dark noise of the SPARC4 cameras.

    """

    def __init__(self, channel: int) -> None:
        """Initialize the class.

        Parameters
        ----------

        channel: integer
            The channel related to the camera.
        """
        self.channel = channel
        self.spreadsheet_path = os.path.join(
            "AIS", "Noise", "spreadsheet", f"Channel {channel}"
        )

    def calculate_read_noise(self, ccd_operation_mode: dict) -> float:
        """Calculate the read noise of the CCD.

        For the conventional mode, it is used the read noise values of the
        Read_noise_and_gain_values spreadsheet

        For the EM mode, the read noise is obtained through an interpolation of
        the values presente by the respective spreadsheet, as a function of the
        CCD EM gain.

        Parameters
        ----------
        ccd_operation_mode: dictionary
            A dictionary with the parameters of the CCD operation mode.

            em_mode : ['EM', 'Conv']
                Electron Multiplying mode of the camera.
            em_gain : float
                EM gain of the camera.
            readout : [0.1, 1, 10, 20, 30]
                Readout rate of the pixels in MHz.
            preamp : [1, 2]
                Pre-amplifier gain.
            binn : [1, 2]
                Binning of the pixels.
        """
        self.em_mode = ccd_operation_mode["em_mode"]
        self.em_gain = ccd_operation_mode["em_gain"]
        self.readout = ccd_operation_mode["readout"]
        self.preamp = ccd_operation_mode["preamp"]
        self.binn = ccd_operation_mode["binn"]

        if self.em_mode == "Conv":
            self._calculate_read_noise_conventional_mode()
        if self.em_mode == "EM":
            self._calculate_read_noise_em_mode()
        return self.read_noise

    def _calculate_read_noise_conventional_mode(self) -> None:
        idx_tab = 1
        if self.readout == 0.1:
            idx_tab += 4
        if self.preamp == 2:
            idx_tab += 2
        idx_tab += self.binn - 1

        ss_name = os.path.join(
            self.spreadsheet_path,
            "read_noise.csv",
        )
        spreadsheet = pd.read_csv(ss_name)
        self.read_noise = float(spreadsheet["Noise"][idx_tab])

    def _calculate_read_noise_em_mode(self) -> None:
        tab_name = os.path.join(
            self.spreadsheet_path,
            "RN_PA"
            + str(int(self.preamp))
            + "B"
            + str(int(self.binn))
            + "HSS"
            + str(int(self.readout))
            + ".csv",
        )
        spreadsheet = pd.read_csv(tab_name, dtype=np.float64)
        column_em_gain = spreadsheet["EM Gain"]
        column_noise = spreadsheet["Noise (e-)"]
        f = interp1d(column_em_gain, column_noise)
        read_noise = f(self.em_gain)

        self.read_noise = float(read_noise)

    def calculate_dark_current(self, temp: float) -> float:
        """Calculate the dark noise.

        Parameters
        ----------
        temp: float
            Temperature, in celsius degree, of the CCD.
            The provided temperature must be in the intevalo of -70 ºC up to -30 ºC.

        Returns
        -------
        dark_current: float
            Dark current, in e-/pix/s, for the respective SPARC4 channel;
        """
        if not (-70 <= temp <= -30):
            raise ValueError(
                f"Expected value for the CCD temperature is outside the range [-70, -30]: {temp}."
            )

        if self.channel == 1:
            dark_current = 24.66 * np.exp(0.0015 * temp**2 + 0.29 * temp)
        elif self.channel == 2:
            dark_current = 35.26 * np.exp(0.0019 * temp**2 + 0.31 * temp)
        elif self.channel == 3:
            dark_current = 9.67 * np.exp(0.0012 * temp**2 + 0.25 * temp)
        elif self.channel == 4:
            dark_current = 5.92 * np.exp(0.0005 * temp**2 + 0.18 * temp)
        else:
            pass
        return dark_current
