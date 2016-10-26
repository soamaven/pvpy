import numpy as np
from scipy import constants, interpolate
from scipy.optimize import fmin
from PowerSpectrum import PowerSpectrum, PhotocurrentSpectrum
import warnings


def ev_to_nm(ev):
    return constants.h * constants.c / (ev * constants.physical_constants["electron volt"][0]) * 1e9


class SolarCell(object):
    def __init__(self, bandgap=1.1, tilt=0, degrees=True, back_reflector=True, celltemp=300):
        """

        :param bandgap:
        :param tilt:
        :param degrees:
        :param celltemp:
        :param back_reflector:
        """
        if bandgap < 0.2:
            warnings.warn("These results will be inaccurate for bandgaps smaller than 0.2eV.")
        self.bandgap = bandgap
        self.bandgap_lambda = ev_to_nm(bandgap)
        if degrees:
            self.tilt = tilt * constants.pi / 180
        else:
            self.tilt = tilt
        self.celltemp = celltemp
        self.ideality = 1
        self.Vc = constants.k * self.celltemp * self.ideality / constants.e
        start_w = 280
        stop_w = int(np.around(self.bandgap_lambda))
        if back_reflector:
            self.solid_emission_angle = PowerSpectrum.solid_angle(180, 180)
        if not back_reflector:
            self.solid_emission_angle = PowerSpectrum.solid_angle(180, 360)
        # Default illumination spectrum is ambient blackbody
        self.illuminationspectrumAmbient = PhotocurrentSpectrum(start_w=start_w, stop_w=stop_w,
                                                                spectra="BlackBody",
                                                                bbtemp=self.celltemp,
                                                                solidangle=self.solid_emission_angle)
        # give the cell an initial perfect absorptivity for SQ type analysis
        self.absorptivity = self.illuminationspectrumAmbient.get_spectrum()
        self.absorptivity[:, 1] = np.ones(self.absorptivity[:, 1].shape)
        # This is an artifact of the ASTM spectrums having nonuniform bins, we have to set
        # the cell band gap to its nearest wavelength
        self.bandgap_lambda = self.absorptivity[-1, 0]
        # give the cell a black body spectrum
        self.cell_bb_spectrum = PhotocurrentSpectrum(start_w=start_w, stop_w=stop_w, bbtemp=self.celltemp,
                                                     spectra="BlackBody",
                                                     solidangle=self.solid_emission_angle)
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.voltage = 0
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        # Initilize the IQE as perfect by copying the unity absorptivity (i.e. normalized absorption)
        self.IQE = self.absorptivity.copy()
        self.LED_eff = 1
        self.generation = 0
        self.illuminationspectrum = self.illuminationspectrumAmbient
        self.set_illumination(self.illuminationspectrumAmbient)
        self.incident_power = self.illuminationspectrum.get_incident_power()
        self.j_nonrad = 0  # Units of amp/m^2
        self.Voc = self.get_Voc()
        self.Isc = self.get_current()
        # Initilize with some empirical guesses
        # use initial empirical guess of Martin Green
        self.ff = (self.Voc/self.Vc - np.log(self.Voc/self.Vc + .72)) / (self.Voc/self.Vc + 1)
        self.Vmpp = self.ff * self.Isc
        self.Impp = self.ff * self.Voc
        self.maxpower = self.Vmpp * self.Impp

    def set_absorptivity(self, alpha):
        self.absorptivity = alpha
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        return

    def __set_generation(self, photocurrentspec):
        # weight the illumination spectrum by the cells probability of absorbing light, and the probability that
        # it will generate an electron-hole pair
        weight = self.absorptivity.copy()
        weight[:, 1] = self.absorptivity[:, 1] * self.IQE[:, 1] * np.cos(self.tilt)
        weightfun = interpolate.interp1d(weight[:, 0], weight[:, 1], kind='linear')
        weight = photocurrentspec.sub_spectrum(weight[0, 0], weight[-1, 0])
        weight[:, 1] = weightfun(weight[:, 0])
        photocurrentspec.weight_spectrum(weight)
        self.generation = photocurrentspec.integrate()
        return self.generation

    def luminescence_spectrum(self, v=0):
        # TODO: Determine if the chemical potential dependent blackbody spectrum will work for non-ideal diodes
        luminescence_spectrum = PhotocurrentSpectrum(start_w=280, stop_w=self.bandgap_lambda,
                                                     bbtemp=self.celltemp,
                                                     spectra="BlackBody",
                                                     solidangle=self.solid_emission_angle,
                                                     v=v)
        luminescence_spectrum.weight_spectrum(self.absorptivity)
        return luminescence_spectrum

    def get_saturation_current(self):
        return self.cell_bb_spectrum.integrate()

    def set_voltage(self, v=0):
        self.luminescencespectrum = self.luminescence_spectrum(v)
        self.voltage = v
        return

    def set_illumination(self, illuminationspectrum):
        self.illuminationspectrum = illuminationspectrum
        self.__set_generation(self.illuminationspectrum)
        # TODO: Figure out when to set incident power. If user puts in a spectrum that is only defined to the bandgap
        # then the "incident power" is lower than it should be. This is avoided by just using the default start_w stop_w
        # for incident spectrums, but perhaps should require the user explicitly input?
        # self.incident_power = illuminationspectrum.get_incident_power()
        self.Voc = self.get_Voc()
        self.Isc = self.get_current()
        return

    def get_current(self, v=None):
        # Following Eqns 3.10, 3.13, and 3.16 of Shockley-Queisser 1961
        # First, get the generation from the illumination spectrum
        Fs = self.generation  # Fs
        # Next, get the blackbody radiation from the cell at 0 voltage
        # aka the non-radiative recombination at v = 0
        Fc0 = self.cell_bb_spectrum  # Fc0
        Fc0 = Fc0.integrate()
        # Below is only valid for ideal rectifier/diode equation
        R0 = Fc0 * (1 - self.LED_eff) / self.LED_eff  # R(0)
        # Next, get the blackbody radiation from the cell at v voltage
        FcV = self.luminescence_spectrum(self.voltage)
        FcV = FcV.integrate()  # FcV
        # Next, get the radiative recombination at v
        Vc = constants.k * self.celltemp / constants.e
        n = self.ideality
        if v is None:
            v = self.voltage
        assert isinstance(v, (int, float, list, np.ndarray)), "Voltage %g is not an int or float." % v
        # TODO: Implement Double Diode Model Here?
        # Rv = R0 * np.exp(v / (Vc * n))  # R(V)
        # print("Fs: %g, FcV: %g, R(0): %g, R(V): %g" % (Fs, FcV, R0, Rv))  # For Debugging
        current = Fs - Fc0 + (Fc0 / self.LED_eff) * (1 - np.exp(v / (Vc * n)))
        return current

    def get_Voc(self):
        external_generation = self.generation
        cell_recombination0 = self.cell_bb_spectrum.integrate()
        Vc = self.Vc
        Voc = Vc * np.log((self.LED_eff * external_generation / cell_recombination0) - self.LED_eff + 1)
        self.Voc = Voc
        return Voc

    def get_max_power(self):
        normvoc = self.Voc / self.Vc
        fffunc = lambda x: -1 * x * self.get_current(v=x) / (self.Voc * self.Isc)
        self.Vmpp = fmin(fffunc, self.Voc * self.ff, disp=False)[0]
        self.Impp = self.get_current(v=self.Vmpp)
        return self.Vmpp * self.Impp

    def get_fill_factor(self):
        self.maxpower = self.get_max_power()
        self.ff = self.maxpower / (self.Voc * self.Isc)
        return self.ff

    def get_effciency(self):
        self.get_fill_factor()
        self.effciency = self.ff * self.Isc * self.Voc / self.incident_power
        return self.effciency
