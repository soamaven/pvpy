# coding=utf-8
import numpy as np

# define constants
q = 1.60217662e-19  # (C) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
Wien = 2.898e-3
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
h = 6.62607004e-34  # (J.s)
c = 299792458.0  # (m s^-1)
hc_q = h * c / q
pi = np.pi

toRad = pi / 180  # use to convert to radians
toDeg = 180 / pi  # use to convert to degrees


# python helpers that are not pv


# basics

def photon_nm2eV(x):
    """ Given wavelength of a photon in um return the energy in eV """
    return hc_q * 1e9 / x


def photon_wavelength(x, units='eV'):
    """ return the energy of photon in eV and joules given teh wave.
    Default return units are in electron-volts.
    """
    energy = h * c / x  # energy in eV, default return value
    if units.upper().lower() == 'eV'.upper().lower():
        return energy
    elif units.upper().lower() is 'joules'.upper().lower():
        energy *= q  # convert to joules
    return energy


def fermi_function(E, Ef, T):
    """ Given the energies in electron volts return the fermi dirac function """
    kt = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kt) + 1.0)
