"""

Channel Creator
===============

This class creates the several channels of the SPARC4. The AIS is based on the
Factory Method to implement the correct star flux calculation as a function of
the optical path of each channel. This package has as abstract class that will
be provided to the Point_Spread_Function Class. Then, the Point_Spread_Function
class will call the correct concrete channel creator for the respective desired
image.
"""

# from sys import exit

import numpy as np

from ..Read_Noise_Calculation import Read_Noise_Calculation
from ..SPARC4_Spectral_Response import (
    Abstract_SPARC4_Spectral_Response,
    Concrete_SPARC4_Spectral_Response_1,
    Concrete_SPARC4_Spectral_Response_2,
    Concrete_SPARC4_Spectral_Response_3,
    Concrete_SPARC4_Spectral_Response_4,
)


class Abstract_Channel_Creator:
    """
    Abstract Channel Creator Class.

    This is an abstracty class that represents the SPARC4 channels. It is
    responsible to calculate the star or the sky flux as a function of the
    properly instrumental response of the SPARC4 channel.

    Parameters
    ----------
    sparc4_operation_mode : [phot, pol]
        Operation mode of the SPARC4: photometric or polarimetric
    """

    def __init__(self, sparc4_operation_mode, wavelength_interval):
        """Initialize the class."""
        self.sparc4_operation_mode = sparc4_operation_mode
        self.wavelength_interval = wavelength_interval
        self._S4_SR = Abstract_SPARC4_Spectral_Response(wavelength_interval)
        self._CHANNEL_ID = 0
        self._SERIAL_NUMBER = 0

    def get_channel_id(self):
        """Return the Channel ID."""
        return f"Channel {self._CHANNEL_ID}"

    def get_serial_number(self):
        """
        Return CCD serial number.

        Returns
        -------
        _SERIAL_NUMBER: 9914, 9915, 9916, 9917
            Serial number of the CCD

        """
        return self._SERIAL_NUMBER

    def calculate_dark_current(self, ccd_temp):
        """
        Calculate the cark current.

        This function calculates the dark current for each SPARC4 CCD. This is
        an abstract function that is extended in the child classes.

        Returns
        -------
        dakr current: float
            Dark current of the respective SPARC4 CCD.

        """
        dark_current = 0
        return dark_current

    def calculate_read_noise(self, ccd_operation_mode):
        """
        Calculate the read noise the CCD.

        The calculation is performed by providing the CCD operation mode to
        the ReadNoiseCalc package
        """
        RN = Read_Noise_Calculation(ccd_operation_mode, self._CHANNEL_ID)
        self.read_noise = RN.calculate_read_noise()

        return self.read_noise

    def apply_sparc4_spectral_response(self, specific_photons_per_second):
        """
        Apply the sparc4 spectral response.

        This function applies the SPARC4 spectral response on the
        calculated star specific flux

        Parameters
        ----------
        star_specific_photons_per_second : array like
            Specific flux of the star

        Returns
        -------

        star_specific_photons_per_second : array like
            The output photons per second per wavelength after applying the spectral response of the
            SPARC4 instrument. If the used SPARC4 operation mode was 'polarimetric', this variable
            will contain the ordinary and extraordinary rays.

        """
        s4_op_mode = self.sparc4_operation_mode
        self._S4_SR.write_specific_photons_per_second(specific_photons_per_second)
        if s4_op_mode["acquisition_mode"] == "polarimetric":
            self._S4_SR.apply_calibration_wheel(s4_op_mode["calibration_wheel"])
            self._S4_SR.apply_retarder(s4_op_mode["retarder"])
            self._S4_SR.apply_analyser()
        self._S4_SR.apply_collimator()
        self._S4_SR.apply_dichroic()
        self._S4_SR.apply_camera()
        self._S4_SR.apply_ccd()
        self.specific_photons_per_second = (
            self._S4_SR.read_specific_photons_per_second()
        )

        return self.specific_photons_per_second


class Concrete_Channel_1(Abstract_Channel_Creator):
    """
    Concreat Channel Creator Class 1.

    This class calculates the star and/or the sky flux as a function of the
    instrumental response of the SPARC4 Channel 1.
    """

    def __init__(self, sparc4_operation_mode, wavelength_interval):
        """Initialize the Channel 1 Class."""
        self.sparc4_operation_mode = sparc4_operation_mode
        self._CHANNEL_ID = 1
        self._SERIAL_NUMBER = 9914
        self._S4_SR = Concrete_SPARC4_Spectral_Response_1(wavelength_interval)

    def calculate_dark_current(self, ccd_temp):
        """Calculate the dark current.

        This function extends the function of the Abstract Channel Creator
        class. Here, the dark current in e-/AUD for the Channel 1 is
        calculated.


        Returns
        -------
        dark_ current: float
            Dark current in e-/ADU for the CCD of the Channel 1 of the SPARC4
            instrument.
        """
        T = ccd_temp
        self.dark_current = 24.66 * np.exp(0.0015 * T**2 + 0.29 * T)
        return self.dark_current


class Concrete_Channel_2(Abstract_Channel_Creator):
    """
    Concreat Channel Creator Class 2.

    This class calculates the star and/or the sky flux as a function of the
    instrumental response of the SPARC4 Channel 2.
    """

    def __init__(self, sparc4_operation_mode, wavelength_interval):
        """Initialize the Channel 2 Class."""
        self.sparc4_operation_mode = sparc4_operation_mode
        self._CHANNEL_ID = 2
        self._SERIAL_NUMBER = 9915
        self._S4_SR = Concrete_SPARC4_Spectral_Response_2(wavelength_interval)

    def calculate_dark_current(self, ccd_temp):
        """Calculate the dark current.

        This function extends the function of the Abstract Channel Creator
        class. Here, the dark current in e-/AUD for the Channel 2 is
        calculated.


        Returns
        -------
        dark_ current: float
            Dark current in e-/ADU for the CCD of the Channel 1 of the SPARC4
            instrument.
        """
        T = ccd_temp
        self.dark_current = 35.26 * np.exp(0.0019 * T**2 + 0.31 * T)
        return self.dark_current


class Concrete_Channel_3(Abstract_Channel_Creator):
    """
    Concreat Channel Creator Class 3.

    This class calculates the star and/or the sky flux as a function of the
    instrumental response of the SPARC4 Channel 3.
    """

    def __init__(self, sparc4_operation_mode, wavelength_interval):
        """Initialize the Channel 3 Class."""
        self.sparc4_operation_mode = sparc4_operation_mode
        self._CHANNEL_ID = 3
        self._SERIAL_NUMBER = 9916
        self._S4_SR = Concrete_SPARC4_Spectral_Response_3(wavelength_interval)

    def calculate_dark_current(self, ccd_temp):
        """Calculate the dark current.

        This function extends the function of the Abstract Channel Creator
        class. Here, the dark current in e-/AUD for the Channel 3 is
        calculated.


        Returns
        -------
        dark_ current: float
            Dark current in e-/ADU for the CCD of the Channel 1 of the SPARC4
            instrument.
        """
        T = ccd_temp
        self.dark_current = 9.67 * np.exp(0.0012 * T**2 + 0.25 * T)
        return self.dark_current


class Concrete_Channel_4(Abstract_Channel_Creator):
    """
    Concreat Channel Creator Class 4.

    This class calculates the star and/or the sky flux as a function of the
    instrumental response of the SPARC4 Channel 4.
    """

    def __init__(self, sparc4_operation_mode, wavelength_interval):
        """Initialize the Channel 4 Class."""
        self.sparc4_operation_mode = sparc4_operation_mode
        self._CHANNEL_ID = 4
        self._SERIAL_NUMBER = 9917
        self._S4_SR = Concrete_SPARC4_Spectral_Response_4(wavelength_interval)

    def calculate_dark_current(self, ccd_temp):
        """Calculate the dark current.

        This function extends the function of the Abstract Channel Creator
        class. Here, the dark current in e-/AUD for the Channel 4 is
        calculated.


        Returns
        -------
        dark_ current: float
            Dark current in e-/ADU for the CCD of the Channel 1 of the SPARC4
            instrument.
        """
        T = ccd_temp
        self.dark_current = 5.92 * np.exp(0.0005 * T**2 + 0.18 * T)
        return self.dark_current
