# coding=utf-8
from pvpy.pveducation.conversions import *

# Solar radiation


def am_intensity(airmass):
    """
    Calculates the solar spectrum intensity for different air-masses. 
    :param airmass: airmass to use
    :type airmass: float
    :return: direct radiation intensity, global radiation intensity (W/m**2)
    :rtype: tuple
    
    >>> '{0[0]:.5f}, {0[1]:.5f}'.format(am_intensity(1.5))
    '0.94710, 1.04181'
    """
    It = 1.353 * (0.7 ** (airmass ** 0.678))
    Ig = It * 1.1
    return It, Ig


def air_mass(angle):
    """
    Find the air mass for a source (i.e. the sun) at an angle from the earth's surface normal vector.
    :param angle: angle (degrees)
    :type : float
    :return: Air Mass (a. u.)
    :rtype: list
    
    >>> np.round(1 / np.cos(np.radians([60])), decimals=10)
    2
    """
    am = 1 / np.cos(np.radians([angle]))
    return am


def am_shadow(s, h):
    """
    Calculates the air mass using the shadow length, s, of a pole of height, h. 
    :param s: shadow length, also the horizontal edge of a right triangle
    :type s: float_or_int
    :param h: height of a pole casting the shadow of length s
    :type h: float_or_int
    :return: the approximate air mass
    :rtype: float_or_int
    
    >>> '{0:.5f}'.format(am_shadow(np.sqrt(3), 1))
    '2.00000'
    """
    am = np.sqrt(1+(s/h)**2)
    return am


def etr_earth(day_no):
    """retrun extrateresstrial radiation at earth (W/m**2)
    given on the day of the year (day) 
    :param day_no: 
    :type day_no: 
    :return: 
    :rtype: 
    
    
    
    """
    return (1 + .033 * (np.cos(np.radians([(day_no - 2.0) * (360 / 365)])))) * 1353


def H0(x):
    """ return the radiant power density (XXX)
    given distance (m) from sun """
    return 2.8942e25 / (x ** 2)


def blackbody_peak(T):
    """Return the peak wavelenth (nm)
    given the temperature (K)"""
    return 1e9 * wien / T


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
    return 23.45 * np.sin(np.radians([(d - 81) * (360 / 365)]))


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
    azi = np.where(HRA > 0, [2 * pi - azi], azi)
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
