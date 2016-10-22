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
        self.area = area
        if not degrees:
            self.tilt = tilt * 180/constants.pi
        else:
            self.tilt = tilt
        if back_reflector:
            self.solid_emission_angle = PowerSpectrum.solid_angle(180, 180)
        if not back_reflector:
            self.solid_emission_angle = PowerSpectrum.solid_angle(180, 360)
        self.celltemp = celltemp
        # give the cell an initial perfect absorptivity for SQ type analysis
        start_w = 280
        stop_w = int(np.around(self.bandgap_lambda))
        wavelengths = np.arange(start_w, stop_w + 1, dtype=int)
        self.absorptivity = np.vstack((wavelengths, np.ones(wavelengths.shape))).T
        # give the cell a black body spectrum
        self.cell_bb_spectrum = PhotocurrentSpectrum(start_w=start_w, stop_w=stop_w, bbtemp=self.celltemp, spectra="BlackBody",
                                              solidangle=self.solid_emission_angle)
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.voltage = 0
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        # Initilize the IQE as perfect by copying the unity absorptivity (i.e. normalized absorption)
        self.IQE = self.absorptivity.copy()
        self.generation = 0
        # Default illumination spectrum is unconcentrated sunlight
        self.illuminationspectrum = PhotocurrentSpectrum(start_w=start_w, stop_w=stop_w, spectra="Dark")
        self.j_nonrad = 0  # Units of amp/m^2
        self.ideality = 1

    def set_absorptivity(self, alpha):
        self.absorptivity = alpha
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        return

    def set_generation(self, photonspectrum):
        # weight the illumination spectrum by the cells probability of absorbing light, and the probability that
        # it will generate an electron-hole pair
        weight = self.absorptivity.copy()
        weight[:, 1] = self.absorptivity[:, 1] * self.IQE[:, 1] * self.area * np.cos(self.tilt)
        photonspectrum.weight_spectrum(weight)
        self.generation = photonspectrum.integrate()
        return self.generation

    def luminescence_spectrum(self, v=0):
        luminescence_spectrum = self.luminescencespectrum.copy()
        Vc = constants.k * self.celltemp / constants.e
        vweight = luminescence_spectrum.get_spectrum()
        vweight[:, 1] *= np.exp(v / Vc)
        luminescence_spectrum.weight_spectrum(vweight)
        return luminescence_spectrum

    def set_voltage(self, v=0):
        self.luminescencespectrum = self.luminescence_spectrum(v)
        self.voltage = 0
        return

    def set_illumination(self, illuminationspectrum):
        self.illuminationspectrum = illuminationspectrum
        return

    def get_current(self):
        # Following Eqn 3.10 of Shockley-Queisser 1961
        # First, get the generation from the illumination spectrum
        self.set_generation(self.illuminationspectrum)
        external_generation = self.generation
        # Next, get the blackbody radiation from the cell at 0 voltage
        # cell_generation0 = self.luminescence_spectrum()
        # cell_generation0 = cell_generation0.integrate()
        # Next, get the blackbody generation from the cell at v voltage
        cell_generationV = self.luminescence_spectrum(self.voltage)
        cell_generationV = cell_generationV.integrate()
        # Next, get the non-radiative recombination at v = 0
        cell_nonrad_recombination = self.j_nonrad
        # Next, get the radiative recombination at v
        Vc = constants.k * self.celltemp / constants.e
        n = self.ideality
        v = self.voltage
        # TODO: Implement Double Diode Model Here?
        cell_rad_recombination = cell_nonrad_recombination * np.exp(v / (Vc * n))
        current = external_generation - cell_generationV + cell_nonrad_recombination - cell_rad_recombination
        # current = np.vstack((cell_generationV[:, 0], current))
        return current
