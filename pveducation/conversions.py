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
__factors = np.concatenate( (np.arange(-24, -3, 3, ), np.arange(-3, 3), np.arange(3, 25, 3)))
__prefixes = np.array(
    [
        'yocto',
        'zepto',
        'atto',
        'femto',
        'pico',
        'nano',
        'micro',
        'mili',
        'centi',
        'deci',
        '',
        'deca',
        'hecto',
        'kilo',
        'mega',
        'giga',
        'tera',
        'peta',
        'exa',
        'zetta',
        'yotta',
    ]
)
__units_dictionary = {}
for prefix, factor in zip(__prefixes, __factors[::-1]):
    __units_dictionary[prefix] = factor

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


def si_units(fun, *args, **kwargs):
    """
    This is meant to decorate all the functions that output standard unit values so that different units can be
    requested easily.
    :param fun: function that outputs SI units and needs to be converted to prefixed units
    :type fun: function
    :param args: pass all arguments
    :type args: list
    :param kwargs: pass all keyword arguments
    :type kwargs: dict
    :return: converter function that scales fun
    :rtype: function
    """
    #TODO: allow for multiprefix units such as "milliamps per centimeter"
    #TODO: unit agnostic prefixes. Probably a regex filter to figure out that "milliamps" should use the "milli" prefix
    def converter(*args, **kwargs):
        if 'units' in kwargs.keys():
            units = kwargs['units']
        else:
            units = ''
        return fun(*args, **kwargs) * 10.0 ** __units_dictionary[units]

    return converter
