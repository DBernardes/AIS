import os

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

__all__ = ["Noise"]


class Noise:

    def __init__(self, channel: int) -> None:
        """Initialize the class.

        Parameters
        ----------

        channel: [1, 2, 3, or 4].
            The channel ID of the camera.
        """
        self.channel = channel
        self.spreadsheet_path = os.path.join(
            os.path.dirname(__file__),
            "spreadsheet",
            f"Channel {channel}",
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

        Returns
        -------
        float:
            The read noise of the camera as a function of its operation mode.
        """
        self.em_mode = ccd_operation_mode["em_mode"]
        self.em_gain = ccd_operation_mode["em_gain"]
        self.readout = ccd_operation_mode["readout"]
        self.preamp = ccd_operation_mode["preamp"]
        self.binn = ccd_operation_mode["binn"]

        if self.em_mode == "Conv":
            self._calculate_read_noise_conventional_mode()
        elif self.em_mode == "EM":
            self._calculate_read_noise_em_mode()
        else:
            raise ValueError(f"A wrong value was passed to the EM mode: {self.em_mode}")
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
        float:
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
