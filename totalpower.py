import urllib2
from astropy.io import fits

filename = 'eovsa_1-18GHz_sp_20140726_164508.fts'
f = urllib2.urlopen('http://ovsa.njit.edu/fits/20140726/'+filename)
g = open(filename,'wb')
g.write(f.read())
f.close()
g.close()
hdulist = fits.open(filename)
tsys = hdulist[3].data[0][1]
