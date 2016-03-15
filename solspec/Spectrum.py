import numpy as np
from scipy import interpolate, integrate


class Spectrum:
    def __init__(self, start_w: float = 280.0, stop_w: float = 4000.0, spectra: str = "AM1.5G"):
        """
        Initilizer for Spectrum class. Builds custom spectrum if variables are changed.
        :param start_w: starting wavelength in nanometers
        :param stop_w: end wavelength in nanometers
        :param spectra: the name of the spectrum you want to use. AM1.5G is most popular and is the ASTMG173-03 global tilt
        standard, the AM1.5D is the direct+circumsolar standard,
        the AM0Etr spectrum is a non-standard zero air-mass spectrum -- Please compare to the ASTM E490 standard
        :return:
        """
        # the first column should be the wavelength in nanometers, the second is the tilt power density/nm in
        # W/(m**2 nm) = J s^-1 m^-2 nm^-1 = C V m^-2 nm^-1
        spectras = {"AM0Etr": 1, "AM1.5G": 2, "AM1.5D": 3}
        self.spectrum = np.genfromtxt("ASTMG173.csv", delimiter=",", skip_header=2)[:, [0, spectras[spectra]]]
        # build custom spectrum if necessary
        if start_w != 280.0 or stop_w != 4000.0:
            self.__bounds_check(start_w)
            self.__bounds_check(stop_w)
            start_ind = np.where(start_w < self.spectrum[:, 0])[0][0]
            stop_ind = np.where(self.spectrum[:, 0] < stop_w)[0][0]
            self.spectrum = self.spectrum[start_ind:stop_ind, :]
        # build spectrum interpolator
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])

    def __bounds_check(self, w: float):
        """
        Checks that the given wavelength is within the bounds of the Spectrum
        :param w: (float) wavelength in nanometers
        :returns: none
        """
        lowerb = self.spectrum[:, 0][0]
        upperb = self.spectrum[:, 0][-1]
        # See if the wavelength is out of bounds, throw error
        if not lowerb < w < upperb:
            print("Wavelength %0.1f out of ASTMG AM1.5 bounds" % w)
            if w < lowerb:
                raise IndexError("Please use the lower bound of %0.2f." % lowerb)
            elif w > upperb:
                raise IndexError("Please use the upper bound of %0.2f." % upperb)
        else:
            return

    def irradiance_at_wavelength(self, w: float):
        """
        Interpolates the spectrum to give the solar spectral power irradiance at the given wavelength.
        :param w: float, wavelength of interest in nanometers
        :return: power irradiance in units of  W/(m**2 nm)
        """
        self.__bounds_check(w)
        irradiance = self.interp(w)
        return irradiance

    def power_density(self, start_w: float, stop_w: float):
        """
        Integrates the AM15G solar spectrum to get the power density in a sub-spectrum
        :param start_w: (float) shortest wavelength in nanometers
        :param stop_w: (float) longest wavelength of spectrum in nanometers
        :return power: (float) the inegrated power of the sub-spectrum
        """
        self.__bounds_check(start_w)
        self.__bounds_check(start_w)
        # get the total number of discreet wavelengths
        waves = self.spectrum[:, 0]
        bin_limit = (np.where(waves < stop_w)[0][-1] + 1) - (np.where(start_w < waves)[0][0])
        # integrate the power
        power_dens = integrate.quad(self.interp, start_w, stop_w, full_output=1, limit=bin_limit)[0]
        return power_dens  # Units Watts/meters^2
