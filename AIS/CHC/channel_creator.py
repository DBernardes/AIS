# -*- coding: utf-8 -*-

"""Channel Creator Class.

This class creates the several channels of the SPARC4. The AIS is based on the
Factory Method to implement the correct star flux calculation as a function of
the optical path of each channel. This package has as abstract class that will
be provided to the Point_Spread_Function Class. Then, the Point_Spread_Function
class will call the correct concrete channel creator for the respective desired
image.
"""


from S4_SR import (Abstract_SPARC4_Spectral_Response,
                   Concrete_SPARC4_Spectral_Response_1,
                   Concrete_SPARC4_Spectral_Response_2,
                   Concrete_SPARC4_Spectral_Response_3,
                   Concrete_SPARC4_Spectral_Response_4)


class Abstract_Channel_Creator:
    """Abstract Channel Creator Class."""

    _channel_ID = 0
    _serial_number = 0

    def __init__(self):
        pass

    def _factory_method(self):
        pass

    def get_channel_ID(self):
        """Return the Channel ID."""
        return f"Channel {self._channel_ID}"

    def get_serial_number(self):
        """Return CCD serial number.

        Returns
        -------
        _serial_number: 9914, 9915, 9916, 9917
            Serial number of the CCD

        """
        return self._serial_number

    def calc_star_flux(self):
        """Calcute the star flux in photons/s."""
        pass


class Concrete_Channel_1(Abstract_Channel_Creator):
    """Concreat Channel Creator Class 1."""

    _channel_ID = 1
    _serial_number = 9914

    def _factory_method(self):
        pass


class Concrete_Channel_2(Abstract_Channel_Creator):
    """Concreat Channel Creator Class 2."""

    _channel_ID = 2
    _serial_number = 9915

    def _factory_method(self):
        pass


class Concrete_Channel_3(Abstract_Channel_Creator):
    """Concreat Channel Creator Class 3."""

    _channel_ID = 3
    _serial_number = 9916

    def _factory_method(self):
        pass


class Concrete_Channel_4(Abstract_Channel_Creator):
    """Concreat Channel Creator Class 4."""

    _channel_ID = 4
    _serial_number = 9917

    def _factory_method(self):
        pass
