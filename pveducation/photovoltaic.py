# This Python file uses the following encoding: utf-8
# -*- coding: utf-8 -*-

from .semiconductors import *

import numpy as np

""" Basic photovoltaic functions
Alpha Version 0.02

Requires numpy

Typical solar units are used, NOT SI units.

wavelength (nm)
Energy of  photon (eV)
semiconductor dimensions (cm)
degrees instead of radians.
Temperature of 298.15 K (25 degC) not 300 K

to do is listed as 9999

The first line on all input files is ignored to allow for column headers
# denotes a comment in input files and is ignored.

Contributions by: sgbowden, ?, ? etc
"""

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
toDeg = 180 / pi  # use to convert to degress

# Solar Cells


def implied_v(Δn, N, T=298.15):
    """
    Return voltage given doping and excess carrier concentration.

    :param Δn: excess carrier concentration
    :type Δn: float
    :param N: float
    :type N: doping level
    :param T: Temperature in Kelvin
    :type T: float
    :return: implied voltage of Si solar cell
    :rtype: float
    """
    return Vt(T) * np.log((Δn + N) * Δn / ni_Si(T) ** 2)


def implied_carrier(V, N, T=298.15):
    """
     Given voltage and doping in Si solar cell determine return carrier concentration.
    :param V: voltage
    :type V: float
    :param N: doping level
    :type N: float
    :param T: Temperature in Kelvin
    :type T: float
    :return: excess carrier concentration
    :rtype: float
    """
    Δn = (-N + np.sqrt(N ** 2 + 4 * ni_Si(T) ** 2 * np.exp(V / Vt(T)))) / 2
    return Δn


def j0side(ni, W, N, D, L, S):
    f = (S * np.cosh(W / L) + D / L * np.sinh(W * L)) / (D / L * np.cosh(W * L) + S * np.sinh(W / L))
    return q * ni ** 2 * (f * D / (L * N))


def efficiency(v_oc, i_sc, ff, A=1):
    """Return efficiency
    given voc (volts), Isc in (amps), FF
    also works for Jsc since area of 1 is assumed
    """
    return 1000 * v_oc * i_sc * ff / A


def current2gen(I):
    """return generation (eh pairs/s)
    given current (amps)
    """
    return I / q


# def J0(ni, We, Ne, De, Le, Se, Nb, Wb, Db, Lb, Sb):
#    '''determines J0, the dark saturation current, under the narrow base diode
# condition where L > W.'''
#    Fe = (Se*np.cosh(We/Le)+De/Le*np.sinh(We*Le))/(De/Le*np.cosh(We*Le)+Se*np.sinh(We/Le))
#    Fb = (Sb*np.cosh(Wb/Lb)+Db/Lb*np.sinh(Wb*Lb))/(Db/Lb*np.cosh(Wb*Lb)+Sb*np.sinh(Wb/Lb))
#    J0 = q*ni**2*(Fe*De/(Le*Ne)+ Fb*Db/(Lb*Nb))
#    return J0


def ideal_diode(V, I0, T=298.15):
    """ideal diode equation
    Return the current given the voltage and saturation"""
    return I0 * np.exp(V / Vt(T))


def solarcell_light_current(V, IL, I0, T=298.15, seriesResistance=None, shuntResistance=None, units="amps"):
    """
    Return current (amps) of a solar cell
    given voltage, light generated current, I0
    also works for J0
    :param V: voltage
    :type V: float
    :param IL: light current
    :type IL: float
    :param I0: saturation current
    :type I0: float
    :param T: temperature in kelvin
    :type T: float
    :param seriesResistance: series resistance of a solar cell in the circuit model
    :type seriesResistance: float
    :param shuntResistance: shunt resistance of a solar cell in the circuit model
    :type shuntResistance: float
    :param units: SI units of current output
    :type units: str
    :return: Current in amps
    :rtype: float
    """
    if shuntResistance and seriesResistance:
        current = IL - I0 * np.exp(V / Vt(T)) - V / shuntResistance
    elif shuntResistance:
        current = IL - I0 * np.exp(V / Vt(T)) - V / shuntResistance
    elif seriesResistance:
        #TODO: implement finding current when there is series resistance, remove the error below
        raise NotImplementedError("Series resistance method is not implemented yet.")
    else:
        current = IL - I0 * np.exp(V / Vt(T))
    return current


def I_cell_Rshunt(V, IL, I0, Rshunt, T=298.15):
    """ return current (A) of a solar cell from   """
    Warning("Deprecated. Use solarcell_light_current() instead.")
    return IL - I0 * np.exp(V / Vt(T)) - V / Rshunt


def voc(IL, I0, n=1, T=298.15):
    """
    Return the open circuit voltage, voc, (volts) from IL(A) and I0(A).
    IL and Io must be in the same units, Eg, (A), (mA) etc
    Using (mA/cm**2) uses J0 and JL instead.
    :param IL:
    :type IL:
    :param I0:
    :type I0:
    :param n:
    :type n:
    :param T:
    :type T:
    :return:
    :rtype:
    """
    return n * Vt(T) * np.log(IL / I0 + 1)


def V_cell(I, IL, I0,  T=298.15):
    """return cell voltage (volts)
    given current (amps), IL, I0
    """
    return Vt(T) * np.log((IL - I) / I0 + 1)

def cell_params(V,I):
    # assumes data is sorted on voltage
    Voc = np.interp(0, -I, V)
    Isc = np.interp(0, V, I)
    idx = np.argmax(V*I)
    Vmp = V[idx]
    Imp = I[idx]
    FF = Vmp*Imp/(Voc*Isc)
    return Voc, Isc, FF, Vmp, Imp


'''
def I_cell_lambert(V, Rs=0.01, Rsh=1000, n=1, I0=1e-12, IL=0.035):
    exp_term = Rsh*(Rs*IL+V)/(n*Vt()*(Rs+Rsh))
    exp_mult = (Rs*I0*Rsh)/(n*Vt()*(Rs+Rsh))
    I=-lambertw(exp_mult*np.exp(exp_term))*n*Vt()/Rs + (Rsh*IL-V)/(Rs+Rsh)
    return I.real
'''


def Pf_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp):
    """return the % resistivity loss in a finger
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    return (L ** 2 * Jmp * Sf * resistivity) / (3 * wf * df * Vmp) * 100.0


def Pf_shading(wf, Sf):
    return (wf / Sf) * 100.0


def Pf_sheet(Sf, Jmp, Rsheet, Vmp):
    return (Sf ** 2 * Jmp * Rsheet) / (12 * Vmp) * 100.0


def Pf_total(L, Jmp, Sf, resistivity, Rsheet, wf, df, Vmp):
    """return the % resistivity loss in a finger
    Given:
        L: finger length (cm)
        Jmp: currrent density at the max power point in A/cm2
        Sf: finger spacing (cm)
    """
    Presistivity = Pf_resistivity(L, Jmp, Sf, resistivity, wf, df, Vmp)
    Pshading = Pf_shading(wf, Sf)
    Psheet = Pf_sheet(Sf, Jmp, Rsheet, Vmp)
    return Presistivity + Pshading + Psheet, Presistivity, Pshading, Psheet



def FF(Vmp, Jmp, Voc, Isc):
    return (Vmp * Jmp) / (Voc * Isc)


def FF_ideal(Voc, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        voc - open circuit voltage (volts)
    """
    voc = normalised_Voc(Voc, ideality, T)
    FF0 = (voc - np.log(voc + 0.72)) / (voc + 1)
    return FF0


def normalised_Voc(Voc, ideality, T=298.15):
    return Voc / (ideality * Vt(T))


def FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        voc - open circuit voltage (volts)
    """
    # voc = normalised_Voc(voc, ideality, T)
    RCH = Voc / Isc
    rs = Rseries / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - 1.1 * rs) + (rs ** 2 / 5.4)
    return FF


def FF_Rsh(Voc, Isc, Rshunt, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        voc - open circuit voltage (volts)
    """
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - ((voc + 0.7) * FF0) / (voc * rsh))
    return FF


def FF_RsRsh(Voc, Isc, Rseries, Rshunt, ideality=1, T=298.15):
    voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rsh = Rshunt / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FFRs = FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15)
    FFRsRsh = FFRs * (1 - ((voc + 0.7) * FFRs) / (voc * rsh))
    return FFRsRsh

# silicon material properties

def ni_Si(T=298.15):
    """return the intrinsic carrier concentration of silicon (cm**-3) based on temperature (K) """
    return 9.38e19 * (T / 300) * (T / 300) * np.exp(-6884 / T)


def optical_properties(fname='OpticalPropertiesOfSilicon.txt'):
    """returns an array with the optical properties of silicon
    usage: wavelength, alpha, absDepth, n, k = optical_properties
    """
    wavelength, abs_coeff, nd, kd = np.loadtxt(fname, skiprows=1, unpack=True)
    return wavelength, abs_coeff, nd, kd

def IQE_emitter(ab, Se, Le, De, xj):
    GF = ((Se * Le / De) + ab * Le - np.exp(-ab * xj) * ((Se * Le / De) * np.cosh(xj / Le) + np.sinh(xj / Le))) / (
        (Se * Le / De) * np.sinh(xj / Le) + np.cosh(xj / Le)) - Le * ab * np.exp(-ab * xj)
    QEE = (Le * ab / (ab * ab * Le * Le - 1)) * GF
    return QEE

def IQE_base(ab, xj_Wd, Sr, Wb, Lb, Db):
    """return quantum efficiency of the base
    Given:
    ab - (cm) absorption coefficient
    xj_Wd - junction depth (cm)
    Sr - surface recombination velocity (cm/s)
    Lb - diffusion lenght of minority carrier in the base (cm)
    Db - diffusivity of minority carriers in the base (cm²/Vs)
    """
    GF = (ab * Lb - (
        (Sr * Lb / Db) * (np.cosh(Wb / Lb) - np.exp(-ab * Wb)) + np.sinh(Wb / Lb) + Lb * ab * np.exp(-ab * Wb)) / (
              (Sr * Lb / Db) * np.sinh(Wb / Lb) + np.cosh(Wb / Lb)))
    QEB = (np.exp(-ab * (xj_Wd))) * (Lb * ab / (ab**2 * Lb**2 - 1)) * GF
    return QEB

def IQE(ab, Wd, Se, Le, De, xj, Sr, Wb, Lb, Db):
    GF = (ab * Le)
    GF = ((Se * Le / De) + ab * Le - np.exp(-ab * xj) * ((Se * Le / De) * np.cosh(xj / Le) + np.sinh(xj / Le))) / (
        (Se * Le / De) * np.sinh(xj / Le) + np.cosh(xj / Le)) - Le * ab * np.exp(-ab * xj)
    QEE = (Le * ab / (ab * ab * Le * Le - 1)) * GF
    
    GF = (ab * Lb - (
        (Sr * Lb / Db) * (np.cosh(Wb / Lb) - np.exp(-ab * Wb)) + np.sinh(Wb / Lb) + Lb * ab * np.exp(-ab * Wb)) / (
              (Sr * Lb / Db) * np.sinh(Wb / Lb) + np.cosh(Wb / Lb)))
    QEB = (np.exp(-ab * (xj + Wd))) * (Lb * ab / (ab * ab * Lb * Lb - 1)) * GF
    QED = np.exp(-ab * xj) * (1 - np.exp(-ab * Wd))
    IQEt = QEE + QEB + QED
    return QEE, QEB, QED, IQEt


def convertSRtoQE(wavelength, QE):
    """'converts a QE in units to spectral response
    assumes that the wavelength is in nm"""
    spectral_response = 1239.8 * QE / wavelength
    return spectral_response


def convertQEtoSR(wavelength, spectral_response):
    """convert SR (A/W) to QE (unit 0 to 1)
    assumes that the wavelegth is in  nm"""
    QE = spectral_response * wavelength / 1239.8
    return QE




def lifetime_SRH(N, Nt, Et, σ_n, σ_p, Δn, T=298.15):
    Nv = 31000000000000000000 * (T / 300) ** 1.85
    Nc = 28600000000000000000 * (T / 300) ** 1.58
    Eg = 1.1246
    vth = 11000000 * (T / 300) ** 0.5
    p0 = N
    n0 = (ni_Si(300) ** 2) / N
    τ_n0 = 1 / (Nt * σ_n * vth)
    τ_p0 = 1 / (Nt * σ_p * vth)
    n1 = Nc * np.exp(-Et / Vt())
    p1 = Nv * np.exp((-Et - Eg) / Vt())
    k_ratio = σ_n / σ_p
    τ_SRH = (τ_p0 * (n0 + n1 + Δn) + τ_n0 * (p0 + p1 + Δn)) / (n0 + p0 + Δn)
    return τ_SRH





# write LIV() input v and I or J as x and y. returns *voc, Jsc, FF, MP, bastardized Rsh, bastardized Rs
