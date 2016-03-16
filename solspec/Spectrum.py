import numpy as np
from scipy import interpolate, integrate, constants


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
        self.start_w = start_w
        self.stop_w = stop_w
        # build custom spectrum if necessary
        if start_w != 280.0 or stop_w != 4000.0:
            self.__bounds_check(start_w)
            self.__bounds_check(stop_w)
            start_ind = np.where(start_w <= self.spectrum[:, 0])[0][0]
            stop_ind = np.where(self.spectrum[:, 0] <= stop_w)[0][-1] + 1
            self.spectrum = self.spectrum[start_ind:stop_ind, :]

    def __bounds_check(self, *wavelengths: float):
        """
        Checks that the given wavelength is within the bounds of the Spectrum
        :param wavelengths: (float) wavelength in nanometers
        :returns: none
        """
        lowerb = self.spectrum[:, 0][0]
        upperb = self.spectrum[:, 0][-1]
        # See if the wavelength(s) is out of bounds, throw error
        for w in wavelengths:
            if not lowerb <= w <= upperb:
                print("Wavelength %0.2f nm out of spectra bounds" % w)
                if w < lowerb:
                    raise IndexError("Please use the lower bound of %0.2f nm." % lowerb)
                elif w > upperb:
                    raise IndexError("Please use the upper bound of %0.2f nm." % upperb)
            else:
                pass
        return

    def irradiance_at_wavelength(self, *wavelengths: float):
        """
        Interpolates the spectrum to give the solar spectral power irradiance at the given wavelength(s).
        :param: wavelengths: float or list of floats, wavelength of interest in nanometers
        :returns: list of power irradiance in units of  W/(m**2 nm)
        """
        self.__bounds_check(*wavelengths)
        for w in wavelengths:
            irradiance = float(self.interp(w))
            yield irradiance

    def power_flux_spectrum(self, *w: float):
        """
        Returns the discrete spectral power irradiance spectrum
        :param w: (floats) list of len(w) = 2 where w[0] is the start wavelength and w[1] is the stop wavelength. Disregards
        any more than the first two args
        :return: (ndarray) the discrete power flux spectrum
        """
        # TODO: Consider how to use interpolation to allow for input of arbitrary wavelengths to be queried quickly
        if not w:
            return self.spectrum.copy()
        else:
            assert len(w) >= 2, "Not enough input wavelengths."
            start_w = w[0]
            stop_w = w[1]
            self.__bounds_check(*w[0:1])
        start_ind = np.where(start_w <= self.spectrum[:, 0])[0][0]
        stop_ind = np.where(self.spectrum[:, 0] <= stop_w)[0][-1] + 1
        return self.spectrum[start_ind:stop_ind, :].copy()

    def photon_flux_spectrum(self, *w: float):
        """
        Returns the discrete spectral photon flux spectrum
        :param w: (floats) list of len(w) = 2 where w[0] is the start wavelength and w[1] is the stop wavelength.
        Disregards any more than the first two args
        :return spec: (ndarray) the discrete photon flux spectrum
        """
        # TODO: Consider how to use interpolation to allow for input of arbitrary wavelengths to be queried quickly
        if not w:
            spectrum = self.spectrum.copy()
        else:
            assert len(w) >= 2, "Not enough input wavelengths."
            start_w = w[0]
            stop_w = w[1]
            self.__bounds_check(*w[0:1])
            start_ind = np.where(start_w <= self.spectrum[:, 0])[0][0]
            stop_ind = np.where(self.spectrum[:, 0] <= stop_w)[0][-1] + 1
            spectrum = self.spectrum[start_ind:stop_ind, :].copy()
        # convert from power to photons; the power bin is 1nm in the visible wavelengths (per the loaded ASTM standard)
        spectrum[:, 1] = spectrum[:, 1] * (spectrum[:, 0] * 1e-9 / (constants.c * constants.h))
        return spectrum

    def integrated_power_flux(self, *w: float):
        """
        Integrates the AM15G solar spectrum to get the power density in a sub-spectrum. Default is the full width of
        Spectrum
        :param w: (floats) list of len(w) = 2 where w[0] is the start wavelength and w[1] is the stop wavelength
        :return power_f: (ndarray) the integrated power of the sub-spectrum
        """
        # get a power spectrum to integrate
        spectrum = self.power_flux_spectrum(*w)
        # create power flux interpolator
        interp_power = interpolate.interp1d(spectrum[:, 0],spectrum[:, 1])
        # get the total number of discreet wavelengths as a bin limit
        bin_limit = len(spectrum[:, 0])
        # integrate the power
        power_f = integrate.quad(interp_power, spectrum[0, 0], spectrum[-1, 0], full_output=1, limit=bin_limit)
        return power_f[0]  # Units Watts/meters^2

    def integrated_photon_flux(self, *w: float):
        """
        Integrates the AM15G solar spectrum to get the integrated photon density in the spectrum.
        Default is the full width of Spectrum
        :param w: (floats) list of len(w) > 2 where w[0] is the start wavelength and w[1] is the stop wavelength
        :return photon_f: (ndarray) the integrated power of the sub-spectrum
        """
        # get a spectrum of photons to integrate
        spectrum = self.photon_flux_spectrum(*w)
        # create photon flux interpolator
        interp_photon = interpolate.interp1d(spectrum[:, 0], spectrum[:, 1])
        # get the total number of discreet wavelengths
        bin_limit = len(spectrum[:, 0])
        # integrate the power
        photon_f = integrate.quad(interp_photon, spectrum[0, 0], spectrum[-1, 0], full_output=1, limit=bin_limit)
        return photon_f[0]  # Units photons/(s meters^2)
