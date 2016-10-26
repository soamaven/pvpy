# coding=utf-8
from copy import copy
import numpy as np
from scipy import interpolate, integrate, constants
from os import path

#Solar Solid angular diameter from earth surface is 1919 arcseconds = 0.009303575 radians
# Source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html accessed on on October 19 2016
SOLAR_SOLID_ANGLE = 2 * constants.pi * (1 - np.cos(.009303575/2))  # ~6.798e-05 steradians

class PowerSpectrum(object):
    def __init__(self, start_w=280.0, stop_w=4000.0, spectra="AM1.5G", bbtemp=5800, mediumrefindex=1,
                 solidangle=SOLAR_SOLID_ANGLE, userspectrum=None, v=0):
        """
        Initilizer for PowerSpectrum class. Builds custom spectrum if variables are passed when creating instance.
        :param start_w: shortest wavelength in nanometers
        :param stop_w: longest wavelength in nanometers
        :param spectra: the name of the spectrum you want to use. AM1.5G is most popular and is the ASTMG173-03
        global tilt standard,
        the AM1.5D is the direct+circumsolar standard,
        the AM0Etr spectrum is a non-standard zero air-mass spectrum -- Please compare to the ASTM E490 standard
        :param bbtemp: Temperature of the blackbody, default is about that of the Sun
        :param mediumrefindex: refractive index of the medium surrounding the blackbody
        :param solidangle: solid angle of the black body seen by the viewer. E.g. the sun's solid angle on the earth is
        the default value, ~5.79e-5 steradians. The solid angle seen by two parallel plates is 2pi. A solar cell with
        a source encompassing the whole sphere it can see is 4pi.
        :return:
        """
        self.bbtemp = bbtemp
        self.mediumrefindex = mediumrefindex
        self.solidangle = solidangle
        self.v = v
        self.start_w = start_w
        self.stop_w = stop_w
        # the first column should be the wavelength in nanometers, the second is the tilt power density/nm in
        # W/(m**2 nm) = J s^-1 m^-2 nm^-1 = C V m^-2 nm^-1
        if userspectrum is not None:
            spectra = "User"
        spectras = {
            "AM0Etr": 1,
            "AM1.5G": 2,
            "AM1.5D": 3,
            "BlackBody": 4,
            "User": 5
        }
        spectra_ind = spectras[spectra]
        if spectra_ind in range(4):
            self.spectrum = np.genfromtxt(path.join(path.dirname(__file__), './ASTMG173.csv'), delimiter=",",
                                          skip_header=2)[:, [0, spectra_ind]]
            if start_w != 280.0 or stop_w != 4000.0:
                self.spectrum = self.sub_spectrum(start_w, stop_w)
                print(self.spectrum[-1])
        elif spectra == "BlackBody":
            self.start_w = start_w
            self.stop_w = stop_w + 1
            self.spectrum = self.blackbody_spectrum(mediumrefindex, solidangle, bbtemp)
        elif spectra == "User":
            assert isinstance(userspectrum, np.ndarray), "Weight spectrum is not a 2D numpy array."
            self.spectrum = userspectrum
        # create the PowerSpectrum interpolator
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])

    def blackbody_spectrum(self, mediumrefindex=1, solidangle=4*constants.pi, bbtemp=5400, v=0):
        """
        Creates a blackbody distribution of power flux for a given solid angle. Default is entire sphere for the sun.
        :param mediumrefindex: refractive index of medium surrounding blackbody
        :param solidangle: solid angle that source subtends to viewer's perspective
        :param bbtemp: temperature of the blackbody in Kelvin
        :return: Nx2 numpy array with a black body spectrum for the spectrum object with start_w and stop_w as limits
        """
        self.solidangle = solidangle
        self.mediumrefindex = mediumrefindex
        self.bbtemp = bbtemp
        # Initilize wavelengths
        wavelengths = np.arange(self.start_w, self.stop_w, dtype=int)
        spectrum = np.vstack((wavelengths, np.ones(wavelengths.shape))).T
        wavelengths = np.arange(self.start_w, self.stop_w, dtype=float) * 1e-9
        # 2n^2hc^2/lambda^5*(exp(-hc/k lambda T)-1) gives units of W sr^-1 m^-3
        numerator = 2 * (self.mediumrefindex ** 2) * constants.h * constants.c ** 2
        exponential = np.exp((constants.h * constants.c - v) / (constants.k * wavelengths * bbtemp))
        spectrum[:, 1] = numerator / ((wavelengths ** 5) * (exponential - 1))
        # Use provided solid angle to and 1m=1-9 get Power Flux per nanometer
        spectrum[:, 1] *= 1e-9 * self.solidangle
        return spectrum

    @staticmethod
    def solid_angle(theta, phi, degrees=True):
        if degrees:
            theta *= constants.pi / 180.
            phi *= constants.pi / 180.
        return phi - phi * np.cos(theta)

    @staticmethod
    def circular_solid_angle(half_angle, degrees=True):
        if degrees:
            half_angle *= constants.pi/180.
        return 2*constants.pi * (1 - np.cos(half_angle))

    def sub_spectrum(self, start_w, stop_w):
        """
        Returns a subset of the PowerSpectrum specified by some bounding wavelengths
        :param start_w: shortest wavelength
        :param stop_w: longest wavelength
        :return: subspec (ndarray) the spectrum between start_w and stop_w
        """
        self.__bounds_check(*[start_w, stop_w])
        start_ind = np.where(start_w <= self.spectrum[:, 0])[0][0]
        stop_ind = np.where(self.spectrum[:, 0] <= stop_w)[0][-1] + 1
        subspec = self.spectrum[start_ind:stop_ind, :].copy()
        return subspec

    def __bounds_check(self, *wavelengths):
        """
        Checks that the given wavelength is between the shortest and longest wavelengths of the PowerSpectrum
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
                    raise IndexError("Wavelength %0.2f too short. Please use the lower bound of %0.2f nm." % (w, lowerb))
                elif w > upperb:
                    raise IndexError("Wavelength %0.2f too large. Please use the upper bound of %0.2f nm." % (w, upperb))
            else:
                pass
        return

    def value_at_wavelength(self, wavelengths):
        """
        Interpolates the spectrum to give the value of the spectrum at the given wavelength(s).
        :param: wavelengths (float, list) wavelength(s) of interest in nanometers
        :returns: values
        """
        self.__bounds_check(*wavelengths)
        for w in wavelengths:
            irradiance = float(self.interp(w))
            yield irradiance

    def integrate(self, *w):
        """
        Integrates the solar spectrum. By default the full width of the spectrum is integrated, but inputting 2 floats
        within the PowerSpectrum bounds will give the integration of sub-spectrum.
        :param w: (floats, ints) shortest and longest wavelengths for a sub-spectrum
        :return power_f: (float) the integrated power of the sub-spectrum
        """
        # deal with subspectrums if necessary
        if not w:
            spectrum = self.spectrum
        else:
            assert len(w) >= 2 and w[0] < w[
                1], 'Error: Too few wavelengths or start wavelength is not shorter than the longest wavelength.'
            spectrum = self.sub_spectrum(w[0], w[1])
        power_f = integrate.trapz(spectrum[:, 1], spectrum[:, 0])
        return power_f  # Units Watts/meters^2

    def get_incident_power(self):
        # Does not include sub-bandgap power
        return self.integrate()

    def get_spectrum(self):
        """
        Returns a copy of the spectrum.
        :return: (ndarray) The discrete spectrum with the wavelengths in [:,0] and the values in [:,1]
        """
        return self.spectrum.copy()

    def weight_spectrum(self, spec_in, kind="linear"):
        """
        Weights a spectrum by a normalized spectrum, e.g. absorption, reflection, transmission at wavelengths in nm
        :param kind: (str or int, optional)interpolation method specification in scipy.interpolat.interp1d:
        Specifies the kind of interpolation as a string (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic, ‘cubic’
        where ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline interpolation of first, second or third order) or as
        an integer specifying the order of the spline interpolator to use. Default is ‘linear’.
        :param spec_in: (np.ndarray) a 2-D array with wavelengths in nm in spec_in[:,0] and normalized vaules in
        spec_in[:,1]
        :return: (np.ndarray) a weighted spectrum in the same format as spec_in
        """
        if spec_in.shape[1] != 2:
            try:
                spec_in = spec_in.transpose()
            except:
                pass
        spec_in = np.squeeze(spec_in)
        assert spec_in.shape[1] == 2, "Weight spectrum is not a 2D numpy array."
        #TODO: Catch errors for when there aren't enough points for default kind, i.e. uncommon case of <4 points in weight
        spec_fun = interpolate.interp1d(spec_in[:, 0], spec_in[:, 1], kind=kind)
        if spec_in[0, 0] >= self.start_w or spec_in[-1, 0] <= self.stop_w:
            self.spectrum = self.sub_spectrum(spec_in[0, 0], spec_in[-1, 0])
        spec_wt = self.spectrum
        spec_wt[:, 1] = spec_fun(spec_wt[:, 0]) * spec_wt[:, 1]
        return

    def copy(self):
        return copy(self)

    def to_PhotonSpectrum(self):
        self.spectrum[:, 1] *= (self.spectrum[:, 0] * 1e-9 / (constants.c * constants.h))
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])
        self.__class__ = PhotonSpectrum
        return

    def to_PhotoCurrentSpectrum(self):
        self.to_PhotonSpectrum()
        self.to_PhotoCurrentSpectrum()
        return


class PhotonSpectrum(PowerSpectrum):
    def __init__(self, start_w=280.0, stop_w=4000.0, spectra="AM1.5G", bbtemp=5800, mediumrefindex=1,
                 solidangle=SOLAR_SOLID_ANGLE, userspectrum=None, v=0):
        """
        Initilizer for PowerSpectrum class. Builds custom spectrum if variables are passed when creating instance.
        :param start_w: shortest wavelength in nanometers
        :param stop_w: longest wavelength in nanometers
        :param spectra: the name of the spectrum you want to use. AM1.5G is most popular and is the ASTMG173-03
        global tilt standard,
        the AM1.5D is the direct+circumsolar standard,
        the AM0Etr spectrum is a non-standard zero air-mass spectrum -- Please compare to the ASTM E490 standard
        :param bbtemp: Temperature of the blackbody, default is about that of the Sun
        :param mediumrefindex: refractive index of the medium surrounding the blackbody
        :param solidangle: solid angle of the black body seen by the viewer. E.g. the sun's solid angle on the earth is
        the default value, ~5.79e-5 steradians. The solid angle seen by two parallel plates is 2pi. A solar cell with
        a source encompassing the whole sphere it can see is 4pi.
        :return:
        """
        super(PhotonSpectrum, self).__init__(start_w=start_w,
                                             stop_w=stop_w,
                                             spectra=spectra,
                                             bbtemp=bbtemp,
                                             mediumrefindex=mediumrefindex,
                                             solidangle=solidangle,
                                             userspectrum=userspectrum,
                                             v=v)
        self.spectrum[:, 1] *= (self.spectrum[:, 0] * 1e-9 / (constants.c * constants.h))
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])

    def to_PowerSpectrum(self):
        self.spectrum[:, 1] /= (self.spectrum[:, 0] * 1e-9 / (constants.c * constants.h))
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])
        self.__class__ = PowerSpectrum
        return

    def to_PhotoCurrentSpectrum(self):
        self.spectrum[:, 1] *= constants.e
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])
        self.__class__ = PhotocurrentSpectrum
        return

    def get_incident_power(self):
        spectrum_copy = self.copy()
        spectrum_copy.to_PowerSpectrum()
        return spectrum_copy.get_incident_power()


class PhotocurrentSpectrum(PhotonSpectrum):
    def __init__(self, start_w=280.0, stop_w=4000.0, spectra="AM1.5G", bbtemp= 5800, mediumrefindex=1,
                 solidangle=SOLAR_SOLID_ANGLE, userspectrum=None, v=0):
        """
        Initilizer for PowerSpectrum class. Builds custom spectrum if variables are passed when creating instance.
        :param start_w: shortest wavelength in nanometers
        :param stop_w: longest wavelength in nanometers
        :param spectra: the name of the spectrum you want to use. AM1.5G is most popular and is the ASTMG173-03
        global tilt standard,
        the AM1.5D is the direct+circumsolar standard,
        the AM0Etr spectrum is a non-standard zero air-mass spectrum -- Please compare to the ASTM E490 standard
        :param bbtemp: Temperature of the blackbody, default is about that of the Sun
        :param mediumrefindex: refractive index of the medium surrounding the blackbody
        :param solidangle: solid angle of the black body seen by the viewer. E.g. the sun's solid angle on the earth is
        the default value, ~5.79e-5 steradians. The solid angle seen by two parallel plates is 2pi. A solar cell with
        a source encompassing the whole sphere it can see is 4pi.
        :return:
        """
        super(PhotocurrentSpectrum, self).__init__(start_w=start_w,
                                                   stop_w=stop_w,
                                                   spectra=spectra,
                                                   bbtemp=bbtemp,
                                                   mediumrefindex=mediumrefindex,
                                                   solidangle=solidangle,
                                                   userspectrum=userspectrum,
                                                   v=v)
        self.spectrum[:, 1] *= constants.e
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])

    def to_PhotonSpectrum(self):
        self.spectrum[:, 1] /= constants.e
        self.interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])
        self.__class__ = PhotonSpectrum
        return

    def to_PowerSpectrum(self):
        self.to_PhotonSpectrum()
        self.to_PowerSpectrum()
        return
