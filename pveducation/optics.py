# coding=utf-8
from pveducation.conversions import *


# Optics
def snell(n1, n2, θ1):
    """ return the refracted angle
    given refractive index of incident medium, transmission medium and incident angle"""
    θ2 = np.degrees(np.arcsin(n1 / n2 * np.sin(np.radians(θ1))))
    return θ2


def extinction2abs(kd, wavelength):
    """absorption coefficienct (cm-1) from extinction coefficient (units)
    what is the wavelength in?"""
    return 1e7 * 4 * π * kd / wavelength


def FresReflect(ni, d):
    """ Return the fresnel reflectivity
    Given refractive index and the dot product XXX
    """
    # gives the amount reflected off a dielectric interface ni incident refractive index, d dotproduct}
    # ni = 1/1.5;
    # d = 0.0;
    c = d
    TI = (1 / ni) ** 2 + c ** 2 - 1
    # for
    ##    if (TI <= 0):
    ##        return 1
    ##    else:
    g = np.sqrt(TI)
    return 0.5 * ((g - c) / (g + c)) ** 2 * (1 + ((c * (g + c) - 1) / (c * (g - c) + 1)) ** 2)


def ARCthick(wavelength, n1):
    """
    return something XXX
    given wavelength and the refractive index
    """
    return wavelength / (4 * n1)


def ARC_opt_n(n0, n2):
    return np.sqrt(n0 * n2)


def ARC_refl(wavelength, n0, n1, nSemi, thickness):
    """
    function to calculate the reflectivity from a semiconductor
    Example call:
    ARC(true,[300:10:1200],1,1.1,2,3.5,100,200);
    n0 - ambient in units
    n1 - refractive index of the dielectric layer 1
    nSemi - refractive index of the semicondcutor
    """

    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - nSemi) / (n1 + nSemi)
    θ = (2 * π * n1 * thickness) / wavelength
    Refl = 100 * (r1 * r1 + r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ)) / (
        1 + r1 * r1 * r2 * r2 + 2 * r1 * r2 * np.cos(2 * θ))
    return Refl


def DLARC_refl(wavelength, n0, n1, n2, nSemi, thickness1, thickness2):
    """return reflectivity (units) of a double layer antireflection coating
    n0 - ambient (units)
    n1 - refractive index of the dielectric layer 1 (units)
    n2 - refractive index of the dielectric layer 2 (units)
    nSemi - refractive index of the semicondcutor
    wavelength, thickness1, thickness 2 all in same units (m) or (nm) etc.
    """
    r1 = (n0 - n1) / (n0 + n1)
    r2 = (n1 - n2) / (n1 + n2)
    r3 = (n2 - nSemi) / (n2 + nSemi)
    θ1 = (2 * π * n1 * thickness1) / wavelength
    θ2 = (2 * π * n2 * thickness2) / wavelength

    numerator = r1 * r1 + r2 * r2 + r3 * r3 + r1 * r1 * r2 * r2 * r3 * r3 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))
    denominator = 1 + r1 * r1 * r2 * r2 + r1 * r1 * r3 * r3 + r3 * r3 * r2 * r2 + 2 * r1 * r2 * (1 + r3 * r3) * np.cos(
        2 * θ1) + 2 * r2 * r3 * (1 + r1 * r1) * np.cos(2 * θ2) + 2 * r1 * r3 * np.cos(
        2 * (θ1 + θ2)) + 2 * r1 * r2 * r2 * r3 * np.cos(2 * (θ1 - θ2))

    return numerator / denominator
