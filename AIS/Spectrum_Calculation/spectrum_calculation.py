"""
Spectrum Calculation Class
===========================

This class calculates the spectrum of an astronomical object as a function of its
magnitude
"""

import numpy as np


class Spectrum_Calculation:
    """Spectrum Calculation class.

    This class calculates the star and the sky spectrum based on the object
    magnitude.
    """

    _H = 6.62607004e-34  # m2 kg / s
    _C = 3e8  # m/s
    _K = 1.38064852e-23  # m2 kg s-2 K-1

    def __init__(self, temperature=5700, l_init=350, l_final=1100, l_step=50):
        """Initialize the class.

        Parameters
        ----------

        temperature: float, optional
            Black-body star temperature in Kelvin.

        l_init: int, optional
            Initial wavelength in nanometers.

        l_final: int, optional
            Final wavelength in nanometers.

        l_step: int, optional
            Step for the wavelength interval in nanometers.

        Returns
        -------

        star_specific_flux: array-like
            Star specific flux.
        """

        self.temperature = temperature
        self.l_init = l_init
        self.l_final = l_final
        self.l_step = l_step

    def calculate_star_specific_flux(self):
        """Calculate the star specific flux.

        This function calculates the star specific flux as a function of the
        provided temperatura. To meet that, it uses the Plank Law:

        :math:`B_{\\lambda}(T) = \\frac{2hc^2}{\\lambda^5} \\frac{1}{e^{hc/\\lambda kT} - 1}`

        where :math:`B_{\\lambda}` is the intensity, :math:`h` is the Plank constante, :math:`c`
        is the light speed, :math:`\\lambda` is the wavelength, :math:`k` is the Bolztman constant,
        and :math:`T` is the temperature. The used wavelength interval is given by the initial, final,
        and step parameters provided to the function.

        """
        T = self.temperature
        h = self._H
        c = self._C
        k = self._K
        specific_flux = []
        for Lambda in range(self.l_init, self.l_final, self.l_step):
            Lambda *= 1e-9
            B = (
                2
                * h
                * c ** 2
                / Lambda ** 5
                * 1
                / (np.e ** (h * c / (Lambda * k * T)) - 1)
            )
            specific_flux.append(B)

        temporary = np.asarray(specific_flux)
        self.specific_flux_length = len(specific_flux)
        specific_flux = np.zeros((4, self.specific_flux_length))
        specific_flux[0, :] = temporary

        self.star_specific_flux = specific_flux

        return self.star_specific_flux

    def calculate_sky_specific_flux(self):
        """Calculate sky specific flux.

        This functions calculates the specific flux of the sky.
        This flux correspond to 10 % of the star flux.
        """

        self.sky_specific_flux = self.calculate_star_specific_flux() * 0.1
        return self.sky_specific_flux
