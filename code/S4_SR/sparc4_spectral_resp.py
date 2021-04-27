"""
SPARC4 Spectral Reponse Class
=============================

This is the SPARC4 Spectral Reponse Class for the calculation of the output
flux of an astronomical object as a funtction of the SPARC4 instrumental
response.
"""


class Abstract_SPARC4_Spectral_Response:
    """Abstract class of the SPARC4 spectral response."""

    _CHANNEL_ID = 0

    def __init__(self, spectrum):
        """Initialize the class.

        Parameters
        ----------
        spectrum : array like
            Spectrum of the astronomical object.

        Returns
        -------
        None.

        """
        self.spectrum = spectrum

    def get_channel_ID(self):
        """Return the chanel ID.

        Returns
        -------
        _CHANNEL_ID : [1, 2, 3, 4]
            Channel ID
        """
        return self._CHANNEL_ID

    def calibration_wheel(self):
        """Calibration wheel spectrum response.

        Apllies the calibration wheel spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def retarder(self):
        """Retarder spectrum response.

        Apllies the retarder spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def analyzer(self):
        """Analyzer spectrum response.

        Apllies the analyzer spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def collimator(self):
        """Collimator spectrum response.

        Apllies the collimator spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def dichroic(self):
        """Dichroic spectrum response.

        Apllies the dichroic spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def camera(self):
        """Camera spectrum response.

        Apllies the camera spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def ccd(self):
        """CCD spectrum response.

        Apllies the CCD spectral response on the flux.
        """
        self.spectrum = self.spectrum

    def integrate_spectrum(self):
        """Integrate the object spectrum.

        Returns
        -------
        flux : float
            The integrated object spectrum as a function of the instrumetal
            response.
        """
        flux = sum(self.spectrum)
        return flux


class Concrete_SPARC4_Spectral_Response_1(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 1."""

    _CHANNEL_ID = 1

    pass


class Concrete_SPARC4_Spectral_Response_2(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 2."""

    _CHANNEL_ID = 2

    pass


class Concrete_SPARC4_Spectral_Response_3(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 3."""

    _CHANNEL_ID = 3
    pass


class Concrete_SPARC4_Spectral_Response_4(Abstract_SPARC4_Spectral_Response):
    """Concrete SPARC4 spectral response of the channel 4."""

    _CHANNEL_ID = 4
    pass
