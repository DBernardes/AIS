"""
Spectrum Calculation Class
===========================

This class calculates the spectrum of an astronomical object as a function of its
magnitude
"""

import numpy as np


class Spectrum_Calculation:
    """
    Spectrum Calculation class.

    This class calculates the star and the sky spectrum based on the object
    magnitude.
    """

    _H = 6.62607004e-34  # m2 kg / s
    _C = 3e8  # m/s
    _K = 1.38064852e-23  # m2 kg s-2 K-1
    _TELESCOPE_EFFECTIVE_AREA = 0.804  # m2
    _SOLAR_DISTANCE = 1.4960e11  # m
    _SOLAR_RADIUS = 6.9551e8  # m
    _BAND_PASS = 0.2e-6  # m
    _S_0 = 4e-2  # W/m2/m

    def __init__(
        self,
        wavelength_interval,
        star_temperature=5700,
        star_radius=1,
        star_dist=1,
    ):
        """
        Initialize the class.

        Parameters
        ----------

        wavelength_interval: array like
            Array with the wavelengths of the spectrum of the star.
        temperature: float, optional
            Black-body star temperature in Kelvin.
        star_radius: float, optional
            Radius of the star in solar unitis.
        star_dist: float, optional
            Distance of the star in solar units.

        Returns
        -------

        star_specific_flux: array-like
            Star specific flux.
        """

        self.star_temperature = star_temperature
        self.wavelength_interval = wavelength_interval
        self.star_radius = star_radius * self._SOLAR_RADIUS
        self.star_dist = star_dist * self._SOLAR_DISTANCE

    def calculate_specific_flux(self, magnitude):
        """Calculate the star specific flux."""
        h = self._H
        c = self._C
        S_0 = self._S_0
        B = self._BAND_PASS
        tel_area = self._TELESCOPE_EFFECTIVE_AREA
        temp = []
        for Lambda in self.wavelength_interval:
            Lambda *= 1e-9
            photons_number = (
                S_0 * 10 ** (-magnitude / 2.5) * Lambda * B * tel_area / (h * c)
            )
            temp.append(photons_number)

        specific_flux = np.zeros((4, len(self.wavelength_interval)))
        specific_flux[0, :] = temp

        return specific_flux

    # def calculate_star_specific_flux_1(self):
    #     """Calculate the star specific flux.

    #     This function calculates the star specific flux as a function of the
    #     provided temperatura. To meet that, it uses the Plank Law:

    #     :math:`B_{\\lambda}(T) = \\frac{2hc^2}{\\lambda^5} \\frac{1}{e^{hc/\\lambda kT} - 1}`

    #     where :math:`B_{\\lambda}` is the intensity, :math:`h` is the Plank constante, :math:`c`
    #     is the light speed, :math:`\\lambda` is the wavelength, :math:`k` is the Bolztman constant,
    #     and :math:`T` is the temperature. The used wavelength interval is given by the initial, final,
    #     and step parameters provided to the function.

    #     """
    #     T = self.temperature
    #     h = self._H
    #     c = self._C
    #     k = self._K
    #     specific_flux = []
    #     areas_relation = (
    #         self.star_radius ** 2 * self._TELESCOPE_EFFECTIVE_AREA / self.star_dist ** 2
    #     )
    #     num = (self.l_final - self.l_init) / self.l_step
    #     for Lambda in np.linspace(self.l_init, self.l_final, int(num)):
    #         Lambda *= 1e-9
    #         var1 = 2 * h * c ** 2 / Lambda ** 5
    #         var2 = np.e ** (h * c / (Lambda * k * T)) - 1
    #         black_body = var1 / var2
    #         photon_energy = h * c / Lambda
    #         photons_per_second = black_body * areas_relation / photon_energy
    #         specific_flux.append(photons_per_second * 1e-25)

    #     self.specific_flux_length = len(specific_flux)
    #     self.star_specific_flux = np.zeros((4, self.specific_flux_length))
    #     self.star_specific_flux[0, :] = specific_flux

    #     return self.star_specific_flux

    # def calculate_sky_specific_flux(self):
    #     """Calculate sky specific flux.

    #     This functions calculates the specific flux of the sky.
    #     This flux correspond to 10 % of the star flux.
    #     """

    #     self.sky_specific_flux = self.calculate_star_specific_flux(22) * 0.1
    #     return self.star_specific_flux * 0.1
