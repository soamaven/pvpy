from conversions import *


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