from astropy.io.fits import Header
from astropy.wcs import WCS

h = Header.fromtextfile('header.txt')
wcs = WCS(h)
Xpix = [[1, 1]]
Xfoc = wcs.sip_pix2foc(Xpix, 0)
print wcs.sip_foc2pix(Xfoc, 0)
