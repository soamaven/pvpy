# coding=utf-8
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

'''
User comments
#suggested by Josh nm2ev = photon_nm2eV
'''

# objective is to have as few imports as possible.
import numpy as np


# define constants
q = 1.60217662e-19  # (C) (units go after the comment line)
eV = q
k = 1.38064852e-23  # (J/K)
Wien = 2.898e-3
Stefan_Boltzmann = 5.670367e-08  # (W m^-2 K^-4)
π = np.pi
h = 6.62607004e-34  # (J.s)
c = 299792458.0  # (m s^-1)
hc_q = h * c / q
pi = np.pi

toRad = pi / 180  # use to convert to radians
toDeg = 180 / pi  # use to convert to degress

# python helpers that are not pv


# basics

def photon_nm2eV(x):
    """ Given wavelength of a photon in um return the energy in eV """
    return hc_q * 1e9 / x


def PhotonWavelength(x):
    """ return the energy of photon in eV and joules given teh wave
    """
    energy_eV = h * c / x
    energy_J = energy_eV * q
    return energy_eV, energy_J


def fermi_function(E, Ef, T):
    """ Given the energies in electron volts return the fermi dirac function """
    kT = k * T / q  # in eV
    return 1 / (np.exp((E - Ef) / kT) + 1.0)


# Solar radiation
def am_intensity(airmass):
    """ return radiation intensity (W/m**2)
    given airmass (units) """
    It = 1.353 * (0.7 ** (airmass ** 0.678))
    Ig = It * 1.1
    return It, Ig


def air_mass(angle):
    """ return Air Mass (units)
    given the XXX angle (degrees) """
    return 1 / np.cos(np.radians(angle))


def am_shadow(s, h):
    """ return the air mass (degress) given a shadow length and height """
    am = 9999
    return am


def etr_earth(day_no):
    """retrun extrateresstrial radiation at earth (W/m**2)
    given on the day of the year (day) """
    return (1 + .033 * (np.cos(np.radians((day_no - 2.0) * (360 / 365))))) * 1353


def H0(x):
    """ return the radiant power density (XXX)
    given distance (m) from sun """
    return 2.8942e25 / (x ** 2)


def blackbody_peak(T):
    """Return the peak wavelenth (nm)
    given the temperature (K)"""
    return 1e9 * Wien / T


def blackbody_integrated(T):
    """return Integrated Irradiance (W/m2/str from a blackbody
    given the temperature (K)"""
    return Stefan_Boltzmann * T ** 4


def blackbody_spectrum(T, wavelength):
    """return the spectral intensity at a T (K) and wavelength(um)
    XXX need to figure out the units"""
    F = 3.742 / ((wavelength ** 5) * (np.exp(1.439e4 / (wavelength * T)) - 1))
    return F


def declination(d):
    """return declination anlln ae of sun (degrees) as a function of day """
    return 23.45 * np.sin(np.radians((d - 81) * (360 / 365)))


def equation_of_time(dayNo):
    """return the equation of time (XXX)
    given the day number """
    B = 360.0 / 365.0 * (dayNo - 81.0)
    B = np.radians(B)
    EoT = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.sin(B)
    # print('EoT', EoT)
    return EoT


def time_correction(EoT, longitude, GMTOffset):
    """ minutes """
    LSTM = 15.0 * GMTOffset
    TimeCorrection = 4.0 * (longitude - LSTM) + EoT
    # print('TC',TimeCorrection)
    return TimeCorrection


def elevation(dec, lat, hra):
    """return the elevation angle of the sun
    given declination, latitude and hour angle of sun
    """
    return np.arcsin(np.sin(dec) * np.sin(lat) + np.cos(dec) * np.cos(lat) * np.cos(hra))


def sun_rise_set(LAT, DEC):
    LAT *= toRad
    DEC *= toRad
    A = (np.sin(LAT) * np.sin(DEC)) / (np.cos(LAT) * np.cos(DEC))
    HH = (np.arccos(A)) * toDeg / 15.0
    sunrise = 12.0 - HH - (TimeCorrection / 60.0)
    sunset = 12 + HH - (TimeCorrection / 60.0)
    return sunrise, sunset


def elev_azi(DEC, LAT, LSTM):
    DEC *= toRad
    LAT *= toRad
    HRA = 15.0 * (LSTM - 12.0)
    # print(HRA)
    HRA *= toRad
    elev = np.arcsin((np.sin(DEC) * np.sin(LAT)) + (np.cos(DEC) * np.cos(LAT) * np.cos(HRA)))
    # print('altitude',alt*toDeg)
    azi = np.arccos((np.cos(LAT) * np.sin(DEC) - np.cos(DEC) * np.sin(LAT) * np.cos(HRA)) / np.cos(elev))

    '''
    if (LSTM > 12):
    	azi = 2*pi - azi
    '''
    azi = np.where(HRA > 0, 2 * pi - azi, azi)
    # print('azimuth',azi*toDeg)
    return elev * toDeg, azi * toDeg


def sun_position(dayNo, latitude, longitude, GMTOffset, H, M):
    """ return the positon of the sun as a elevation and azimuth given
    latitude, logitude and the GMTOffset, """
    EoT = equation_of_time(dayNo)
    TimeCorrection = time_correction(EoT, longitude, GMTOffset)
    DEC = declination(dayNo)
    LSTM = H + (TimeCorrection + M) / 60.0
    elev, azi = elev_azi(DEC, latitude, LSTM)
    # print(elev,azi)
    # print(dayNo,DEC,latitude,LSTM, TimeCorrection)
    return elev, azi

#Optics
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



# Basic Semiconductors

def tau_from_L(L, D):
    """
    Return the lifetime (s)
    given the  diffusion length (cm) and diffusivity (cm2/s)
    """
    return np.sqrt(L * D)


def diff_length(lifetime, diffusivity):
    """ return carrier diffusion length (cm)
    given carrier lifetime(s) and diffusivity (units)
    """
    return np.sqrt(lifetime * diffusivity )


def tau_b_from_tau_eff(tau_eff, S, W):
    """Return the bulk lifetime in us
    Given tau_eff (us)
    surface recombination, cm/s
    W, cm
    """
    return 9999


def Vt(T=298.15):
    """return thermal voltage (volts) given temperature(Kelvin).
    Default temp is 298.15 K"""
    return k * T / q

def masetti_mobility(N):
    """ mobility model """
    µmax = 1414
    µmin = 68.5
    u1 = 56.1
    Nref1 = 9.20e16
    Nref2 = 3.41e20
    a = 0.711
    b = 1.98
    return µmin + (µmax - µmin) / (1 + ((N / Nref1) ** a)) - u1 / (1 + ((Nref2 / N) ** b))


def phos_active(T):
    """return the active limit of phosphorous from temperature (degC)"""
    T = sp.constants.C2K(T)
    return 1.3e22 * np.exp(-0.37 * eV / (k * T))


def phos_solubility(T):
    """return the solubility limit of phosphorous from temperature (degC)"""
    T = sp.constants.C2K(T)
    return 2.45e23 * np.exp(-0.62 * eV / (k * T))


def carrier_conc(N, ni=8.6e9):
    """return carrier concentration
    given doping and intrinsic carrier concentration
    Typically (cm) used for all units but (m) etc would also work"""
    majority = N
    minority = N / (ni ** 2)
    return majority, minority


def mob_thurber(N, p_type=True, majority=True):
    """ mobility of electrons as minority carriers in silicon based on doping (cm**-3) """
    i = 2 * p_type + majority
    # n-type minority, n-type majority, p-type minority, p-type majority
    umax = [1417, 1417, 470, 470][i]
    umin = [160, 60, 155, 37.4][i]
    Nref = [5.6e16, 9.64E16, 1e17, 2.82E17][i]
    a = [0.647, 0.664, 0.9, 0.642][i]
    return umin + (umax - umin) / (1 + ((N / Nref) ** a))


def resistivity_Si_n(Ndonor):
    """resistivity of silicon as a function of doping """
    n_minority = ni_Si() ** 2 / Ndonor
    return 1 / ((q * mob_thurber(Ndonor, False) * Ndonor) + (q * mob_thurber(n_minority, False, False) * n_minority))


def resistivity_Si_p(Nacceptor):
    """resistivity of silicon as a function of doping """
    n_minority = ni_Si() ** 2 / Nacceptor
    return 1 / ((q * mob_thurber(Nacceptor) * Nacceptor) + (q * mob_thurber(n_minority, True, False) * n_minority))


# Solar Cells

def impliedV(Δn, N, T=298.15):
    """return voltage
    given doping and excess carrier concentration """
    return Vt(T) * np.log((Δn + N) * Δn / ni_Si(T) ** 2)


def implied_carrier(V, N, T=298.15):
    """ return carrier concentration
    Given voltage and doping determine """
    Δn = (-N + np.sqrt(N ** 2 + 4 * ni_Si(T) ** 2 * np.exp(V / Vt(T)))) / 2
    return Δn


def J0side(ni, W, N, D, L, S):
    F = (S * np.cosh(W / L) + D / L * np.sinh(W * L)) / (De / L * np.cosh(W * L) + S * np.sinh(W / L))
    return q * ni ** 2 * (F * D / (L * N))


def efficiency(Voc, Isc, FF, A=1):
    """Return efficiency
    given Voc (volts), Isc in (amps), FF
    also works for Jsc since area of 1 is assumed
    """
    return 1000 * Voc * Isc * FF / A


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


def I_diode(V, I0, T=298.15):
    """ideal diode equation
    Return the current given the voltage and saturation"""
    return I0 * np.exp(V / Vt(T))


def I_cell(V, IL, I0, T=298.15):
    """ return current (amps) of a solar cell
    given voltage, light generated current, I0
    also works for J0
    """
    return IL - I0 * np.exp(V / Vt(T))


def I_cell_Rshunt(V, IL, I0, Rshunt, T=298.15):
    """ return current (A) of a solar cell from   """
    return IL - I0 * np.exp(V / Vt(T)) - V / Rshunt


def Voc(IL, I0, n=1,  T=298.15):
    """return the open circuit voltage, Voc, (volts) from IL(A) and I0(A).
    IL and Io must be in the same units, Eg, (A), (mA) etc
    Using (mA/cm**2) uses J0 and JL instead.
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
        Voc - open circuit voltage (volts)
    """
    voc = normalised_Voc(Voc, ideality, T)
    FF0 = (voc - np.log(voc + 0.72)) / (voc + 1)
    return FF0


def normalised_Voc(Voc, ideality, T=298.15):
    return Voc / (ideality * Vt(T))


def FF_Rs(Voc, Isc, Rseries, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
    """
    # voc = normalised_Voc(Voc, ideality, T)
    RCH = Voc / Isc
    rs = Rseries / RCH
    FF0 = FF_ideal(Voc, ideality, T)
    FF = FF0 * (1 - 1.1 * rs) + (rs ** 2 / 5.4)
    return FF


def FF_Rsh(Voc, Isc, Rshunt, ideality=1, T=298.15):
    """Return the FF (units)
    Given:
        Voc - open circuit voltage (volts)
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





# write LIV() input v and I or J as x and y. returns *Voc, Jsc, FF, MP, bastardized Rsh, bastardized Rs
