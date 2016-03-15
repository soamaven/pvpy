import numpy as np
from scipy import interpolate

class ASTMGspectrum:
    def __init__(self):
        # the first column should be the wavelength in nanometers, the second is the power density/nm in
        # W/(m**2 nm) = J s^-1 m^-2 nm^-1 = C V m^-2 nm^-1
        self.spectrum = np.genfromtxt("ASTMG173.csv", delimiter=",", skip_header=2)[:, [0, 2]]
        self.AM15interp = interpolate.interp1d(self.spectrum[:, 0], self.spectrum[:, 1])