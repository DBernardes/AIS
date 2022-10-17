"""
Spectral Respose Class
======================


The Spectral Response is an abstract class that represents the optical systems in the light path. 
Theses systems are the atmosphere, the telescope, and the SPARC4 instrument.
"""


import os
import pandas as pd


class Spectral_Response:
    """
    Spectral Response class

    The Spectral Response is an abstract class that represents the optical systems present in the light path 
    between the object and a detector on the ground. 

    Parameters
        ----------

        csv_file_name: str
            The file name containing the spectral response of the optical system.

    Yields
    ------
        spectral_response: array_like
            The spectral response of the optical system.
    """
    _BASE_PATH = os.path.join('AIS', 'Spectral_Response')

    def __init__(self, csv_file_name: str) -> None:
        """Initialize the class."""
        self.csv_file_name = csv_file_name
        return

    @staticmethod
    def _read_csv_file(csv_file_name: str) -> tuple:
        ss = pd.read_csv(csv_file_name)
        return ss['Wavelength (nm)'], ss['Transmitance (%)']

    def get_spectral_response(self) -> list:
        """Return the spectral Response.

        Parameters
        ----------

        spectral_response: array_like.
            The spectral response of the optical system.            
        """
        csv_file_name = os.path.join(self._BASE_PATH, self.csv_file_name)
        _, transmitance = self._read_csv_file(csv_file_name)

        return transmitance
