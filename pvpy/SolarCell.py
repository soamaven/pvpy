from PowerSpectrum import *


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
        stop_w = int(round(self.bandgap_lambda))+1
        self.absorptivity = np.ones((stop_w - start_w, 2))
        self.absorptivity[:, 0] = np.arange(start_w, stop_w, dtype=int)
        # give the cell a black body spectrum
        self.cell_bb_spectrum = PowerSpectrum(start_w, stop_w, bbtemp=self.celltemp, spectra="BlackBody",
                                              solidangle=self.solid_emission_angle)
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        self.voltage = 0
        self.luminescencespectrum = self.cell_bb_spectrum.copy()
        # Initilize the ERE as perfect by copying the unity absorptivity (i.e. normalized absorption)
        self.ERE = self.absorptivity.copy()
        self.generation = 0
        self.illuminationspectrum = None

    def set_absorptivity(self, alpha):
        self.absorptivity = alpha
        self.cell_bb_spectrum.weight_spectrum(self.absorptivity)
        return

    def set_generation(self, photonspectrum):
        self.generation = photonspectrum.integrate() * self.area * np.cos(self.tilt)
        return self.generation

    def luminescence_spectrum(self, v=0):
        luminescence_spectrum = self.luminescencespectrum.copy()
        Vc = constants.k * self.celltemp / constants.e
        vweight = luminescence_spectrum.get_spectrum()
        vweight[:, 1] *= np.exp(v / Vc)
        luminescence_spectrum.weight_spectrum(vweight)
        return luminescence_spectrum.get_spectrum()

    def set_voltage(self, v=0):
        self.luminescencespectrum = self.luminescence_spectrum(v)
        self.voltage = 0
        return

    def set_illumination(self, illuminationspectrum):
        self.illuminationspectrum = illuminationspectrum
        return
