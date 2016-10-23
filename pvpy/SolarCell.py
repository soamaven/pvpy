import numpy as np
from scipy import interpolate, integrate, constants
from PowerSpectrum import PowerSpectrum, PhotocurrentSpectrum, PhotonSpectrum


def ev_to_nm(ev):
    return constants.h * constants.c / (ev * constants.physical_constants["electron volt"][0]) * 1e9


class SolarCell(object):
    def __init__(self, bandgap=1.1, area=1, tilt=0, degrees=True, back_reflector=True, celltemp=300):
        """

        :param bandgap:
        :param area:
        :param tilt:
        :param degrees:
        :param BBtemp:
        :param back_reflector:
        """
        self.bandgap = bandgap
        self.bandgap_lambda = ev_to_nm(bandgap)
        if not degrees:
            self.tilt = tilt * 180/constants.pi
        else:
            self.tilt = tilt
        self.celltemp = celltemp
        start_w = 280
        stop_w = int(np.around(self.bandgap_lambda))
        if back_reflector:
            self.area = area
        if not back_reflector:
            self.area = area * 2
        self.solid_emission_angle = PowerSpectrum.solid_angle(180, 180)
        # give the cell an initial perfect absorptivity for SQ type analysis
        wavelengths = np.arange(start_w, stop_w + 1, dtype=int)
        self.absorptivity = np.vstack((wavelengths, np.ones(wavelengths.shape))).T
        # give the cell a black body spectrum
        self.cell_bb_spectrum = PhotocurrentSpectrum(start_w=start_w, stop_w=20000, bbtemp=self.celltemp,
                                                     spectra="BlackBody",
                                                     solidangle=self.solid_emission_angle)
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        # Default illumination spectrum is ambient blackbody
        self.illuminationspectrumAmbient = PhotocurrentSpectrum(start_w=start_w, stop_w=stop_w,
                                                                spectra="BlackBody",
                                                                bbtemp=self.celltemp,
                                                                solidangle=self.solid_emission_angle)
        self.voltage = 0
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        # Initilize the IQE as perfect by copying the unity absorptivity (i.e. normalized absorption)
        self.IQE = self.absorptivity.copy()
        self.LED_eff = 1
        self.generation = 0

        self.illuminationspectrum = self.illuminationspectrumAmbient
        self.j_nonrad = 0  # Units of amp/m^2
        self.ideality = 1

    def set_absorptivity(self, alpha):
        self.absorptivity = alpha
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        return

    def __set_generation(self, photonspectrum):
        # weight the illumination spectrum by the cells probability of absorbing light, and the probability that
        # it will generate an electron-hole pair
        weight = self.absorptivity.copy()
        weight[:, 1] = self.absorptivity[:, 1] * self.IQE[:, 1] * self.area * np.cos(self.tilt)
        photonspectrum.weight_spectrum(weight)
        self.generation = photonspectrum.integrate()
        return self.generation

    def luminescence_spectrum(self, v=0):
        luminescence_spectrum = self.cell_bb_spectrum.copy()
        Vc = constants.k * self.celltemp / constants.e
        n = self.ideality
        vweight = luminescence_spectrum.get_spectrum()
        vweight[:, 1] *= np.exp(v / (Vc * n))
        luminescence_spectrum.weight_spectrum(vweight)
        return luminescence_spectrum

    def set_voltage(self, v=0):
        self.luminescencespectrum = self.luminescence_spectrum(v)
        self.voltage = v
        return

    def set_illumination(self, illuminationspectrum):
        self.illuminationspectrum = illuminationspectrum
        self.__set_generation(self.illuminationspectrum)
        return

    def get_current(self, v=None):
        # Following Eqns 3.10, 3.13, and 3.16 of Shockley-Queisser 1961
        # First, get the generation from the illumination spectrum
        Fs = self.generation  # Fs
        # Next, get the blackbody radiation from the cell at 0 voltage
        # aka the non-radiative recombination at v = 0
        Fc0 = self.cell_bb_spectrum  # Fc0
        Fc0 = Fc0.integrate()
        R0 = Fc0 * (1 - self.LED_eff) / self.LED_eff  # R(0)
        # Next, get the blackbody radiation from the cell at v voltage
        FcV = self.luminescence_spectrum(self.voltage)
        FcV = FcV.integrate()  # FcV
        # Next, get the radiative recombination at v
        Vc = constants.k * self.celltemp / constants.e
        n = self.ideality
        if v is None:
            v = self.voltage
        assert isinstance(v, (int, float)), "Voltage %g is not an int or float." % v
        # TODO: Implement Double Diode Model Here?
        Rv = R0 * np.exp(v / (Vc * n))  # R(V)
        #print("Fs: %f, FcV: %f, R(0): %f, R(V): %f" % (Fs, FcV, R0, Rv))  # For Debugging
        current = Fs - Fc0 + (Fc0 / self.LED_eff) * (1 - np.exp(v / (Vc * n)))
        return current

    def get_Voc(self):
        self.__set_generation(self.illuminationspectrum)
        external_generation = self.generation
        cell_recombination0 = self.cell_bb_spectrum.integrate()
        Vc = constants.k * self.celltemp / constants.e
        Voc = Vc * np.log((self.LED_eff * external_generation / cell_recombination0) - self.LED_eff + 1)
        return Voc
