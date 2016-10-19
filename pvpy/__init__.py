from os import path
import numpy as np
from pvpy.PowerSpectrum import *
from pvpy import PowerSpectrum

try:
    filename = 'ASTMG173.csv'
    with open(path.join(__path__[0], filename)) as file:
        pass
except IOError as e:
    print("Unable to open file %s" % filename)  # Does not exist OR no read permissions
# TODO: Maybe allow direct download of ASTM1.5G spectrum if it is missing? It shoudln't be missing.
#    import os
#    if not os.path.isfile('ASTMG173.csv'):
#        try:
#            import urllib2, StringIO, tarfile
#            data_url = 'http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/compressed/ASTMG173.csv.tar'
#            download_as_string = urllib2.urlopen(data_url).read()
#            download_as_file = StringIO.StringIO(download_as_string)
#            download_as_tarfile_object = tarfile.open(fileobj=download_as_file)
#            download_as_tarfile_object.extractfile('ASTMG173.csv')
#        except:
#            print("Unable to open file")  # Does not exist OR no read permissions
#            raise


def jsc(*args, **kwargs):
    """
    Gives the absorbed photocurrent in mA/cm**2 of a normalized spectrum.
    :param args: (array like) ideally (N,2)D numpy array with spec_in[:,0] as the wavelengths in nm and spec_in[:,1] as
     the values, but can also be two lists/vectors
    :param kind: (str) interpolation kind. see scipy.interpolate.interp1d
    :param spectra:
    :return:
    """
    kind = kwargs.pop('kind', 'linear')
    spectra = kwargs.pop('spectra', "AM1.5G")
    if len(args) > 1:
        assert args[0].size == args[1].size, \
            "arg1 is %g, arg2 is %g. The inputs should be the same sizes." % (args[0].size, args[1].size)
        args = [np.asarray(arg) if type(arg) is list else arg for arg in args]
        spec_in = np.column_stack((args[0], args[1]))
        assert np.shape(spec_in)[1] == 2, "The input arguments could not be converted to a (N,2) dimensional array."
    elif len(args[0].shape) > 1 and args[0].shape[1] != 1:
        spec_in = args[0].transpose()
    else:
        spec_in = args[0]
    spec_in = np.squeeze(spec_in)
    spec = PowerSpectrum.PhotocurrentSpectrum(spec_in[0, 0], spec_in[-1, 0], spectra)
    spec.weight_spectrum(spec_in, kind=kind)
    return spec.integrate() * .1
