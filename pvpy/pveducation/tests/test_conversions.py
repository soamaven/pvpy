import pytest
from itertools import combinations
from pvpy.pveducation import conversions as cnv
import numpy as np

allunitmods = [k for d in [cnv.__pre_mods.keys(), cnv.__post_mods.keys()] for k in d]
unit_perms = [" ".join(y) for y in combinations(allunitmods, r=3)]
testunits = zip(np.random.rand(unit_perms.__sizeof__()) * 100, unit_perms)


def myfloater(z, *args, **kwargs):
    return z


@pytest.mark.parametrize('x, units',
                         testunits,
                         )
def test_si_units(x, units):
    wrappedfloater = cnv.si_units(myfloater)
    with pytest.raises((KeyError, AssertionError)):
        wrappedfloater(x, units=units)


def gentestunits():
    unit_combs = combinations(cnv.__si_mods.keys(), r=3)
    for comb in unit_combs:
        yield (np.random.rand() * 100, ' '.join(comb))


print("{} tests to do.".format(sum(1 for x in gentestunits())))
unitgen = gentestunits()


@pytest.mark.parametrize('x, units',
                         unitgen,
                         )
def test_si_units2(x, units):
    wrappedfloater = cnv.si_units(myfloater)
    # with pytest.raises((KeyError, AssertionError)):
    wrappedfloater(x, units=units)


def test_nm2ev():
    assert cnv.photon_nm2eV(2.0, units='nanometers') == pytest.approx(619.9209872915753, 1e-16)


def test_photonwavelength():
    assert cnv.photon_wavelength(551.5, e_units='ev', units='nano') == pytest.approx(2.2481268791697064, 1e-16)


def test_fermidirac():
    assert cnv.fermi_function(2.0, 2.02, 300.23) == pytest.approx(0.68417860157406751, 1e-16)
