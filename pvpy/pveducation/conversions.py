# coding=utf-8
import numpy as np
import re
from functools import wraps

# define constants
q = 1.60217662e-19  # (C) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
wien = 2.898e-3
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
h = 6.62607004e-34  # (J.s)
c = 299792458.0  # (m s^-1)
hc_q = h * c / q
pi = np.pi

toRad = pi / 180  # use to convert to radians
toDeg = 180 / pi  # use to convert to degrees
__factors = np.concatenate((np.arange(-24, -3, 3, dtype="float64"),  # up to milli (py2, neg int exponents invalid)
                            np.arange(-3, 0, dtype="float64"),  # up to base (py2, neg int exponents invalid)
                            np.arange(1, 3, dtype="int"),  # base up to hecta
                            np.arange(3, 25, 3, dtype="int")))  # kilo to yotta
__prefixes = np.array(
    [
        'yocto',
        'zepto',
        'atto',
        'femto',
        'pico',
        'nano',
        'micro',
        'milli',
        'centi',
        'deci',
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
__si_mods = {}
for prefix, fact in zip(__prefixes, __factors[::-1]):
    __si_mods[prefix] = fact
# Special english modifiers for units
__pre_mods = {"cubic": 3, "square": 2, "per": -1.0, "inverse": -1.0}
__post_mods = {"cubed": 3, "squared": 2}


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


# decorator to allow user to output preferred units from functions
def si_units(fun):
    """
    This is meant to decorate all the functions that output standard unit values so that different units can be
    requested easily.
    :param fun: function that outputs SI units and needs to be converted to prefixed units
    :type fun: function
    :return: converter function that scales fun
    :rtype: function
    """

    regex_parser = re.compile(r"\w+", )

    @wraps(fun)
    def units_wrapper(*args, **kwargs):
        # check that the caller is trying to modify the units before doing anything
        if "units" in kwargs.keys():
            # remove the "units" keyword so it doesn't get passed to wrapped function
            units = kwargs.pop("units")
            assert isinstance(units, str), "Units should be a sentence like string."
            # Capture all the words in the units string
            matches = regex_parser.findall(units)
            # pre-allocate a list of unity factors
            unit_factors = [1] * len(matches)
            items = []
            # Appropriately apply the factors
            for match in matches:
                # Capture all the factors corresponding to the exponents of 10 needed to change units
                # Will NOT capture odd redundant prefixes or modifiers, e.g. centi-centimeters or per per meter
                for (key, value) in zip(__si_mods.keys(), __si_mods.values()):
                    if key in match:
                        items.append((key, value))
                for (key, value) in zip(__pre_mods.keys(), __pre_mods.values()):
                    if key == match:
                        items.append((key, value))
                        assert match is not matches[0], "{:s} unit modifier can't be first!".format(match)
                for (key, value) in zip(__post_mods.keys(), __post_mods.values()):
                    if key == match:
                        items.append((key, value))
                        assert match is not matches[-1], "{:s} unit modifier can't be last!".format(match)

            for i, (item, match) in enumerate(zip(items, matches)):
                if item[0] in __si_mods.keys():
                    unit_factors[i] *= item[1]
                    continue
                elif item[0] in __pre_mods.keys():
                    try:
                        # Apply the pre modifiers to the next position
                        unit_factors[i + 1] *= (item[1] * unit_factors[i])
                        # If next modifier is also a pre modifier, propagate the modifiers forward
                        if matches[i + 1] in __pre_mods.keys():
                            unit_factors[i + 2] *= unit_factors[i + 1]
                            unit_factors[i + 1] = 1
                        continue
                    except IndexError as e:
                        print(e.message)
                        raise(KeyError, "{:s} are/is invalid unit ordering.".format(" ".join(matches[i:])))
                elif item[0] in __post_mods.keys():
                    try:
                        # Apply the post modifiers to the previous position
                        unit_factors[i - 1] *= (item[1] * unit_factors[i])
                        # If next modifier is also a modifier, propagate the modifiers
                        if matches[i - 1] in __post_mods.keys():
                            unit_factors[i] *= unit_factors[i - 1]
                            unit_factors[i - 1] = 1
                        continue
                    except IndexError as e:
                        print(e.message)
                        raise (KeyError, "{:s} are/is invalid unit ordering.".format(" ".join(matches[0:(i+1)])))
                else:
                    raise (SyntaxError("{:s} can't be used in a units string.".format(match)))
            # Last thing to do is get rid of all the unit factors
            unit_factors[:] = (factor for factor in unit_factors if factor != 1.0)
            exponent = np.sum(unit_factors)
            return fun(*args, **kwargs) * 10 ** exponent
        else:
            return fun(*args, **kwargs)

    return units_wrapper


if __name__ == '__main__':
    # pass
    def myfloater(z):
            return z
    wrappedfloater = si_units(myfloater)
    wrappedfloater(5, units='cubic per square inverse')
