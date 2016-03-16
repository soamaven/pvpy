from solspec import PowerSpectrum
import numpy as np

try:
    with open('ASTMG173.csv') as file:
        pass
except IOError as e:
    print("Unable to open file")  # Does not exist OR no read permissions
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


def jsc(spec_in: np.ndarray, kind: str ="linear", spectra: str ="AM1.5G"):
    """
    Gives the absorbed photocurrent in mA/cm**2 of a normalized spectrum.
    :param spec_in: ((N,2) numpy array) with spec_in[:,0] as the wavelengths in nm and spec_in[:,1] as the values
    :param kind: (str) interpolation kind. see scipy.interpolate.trapz
    :param spectra:
    :return:
    """
    if spec_in.shape[1] != 2:
        try:
            spec_in = spec_in.transpose()
        except:
            pass
    spec_in = np.squeeze(spec_in)
    spec = PowerSpectrum.PhotocurrentSpectrum(spec_in[0, 0], spec_in[-1, 0], spectra)
    spec.weight_spectrum(spec_in, kind=kind)
    return spec.integrate() * .1
