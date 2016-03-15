from Spectrum import Spectrum


#try:
#    with open('ASTMG173.csv') as file:
#        pass
#except IOError as e:
#    print("Unable to open file")  # Does not exist OR no read permissions
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