"""
Testing-out WCS conversions in astropy
"""
import os, sys
sys.path.append(os.path.abspath('..'))

import matplotlib.pyplot as plt
import numpy as np

from astropy.wcs import WCS
from astropy.io import fits

from spheredb.get_data import get_stripe82_file


def construct_wcs_from_file():
    print "Constructing WCS from a FITS file:"
    # Construct a wcs instance from file
    hdulist = get_stripe82_file(301, 2700, 1, 'g', 18)
    wcs = WCS(hdulist[0].header)
    
    # Test a round-trip from pix to world and back
    origin = 0
    pix_coords = np.arange(1, 11).reshape((5, 2))
    world_coords = wcs.wcs_pix2world(pix_coords, origin)
    pix_coords2 = wcs.wcs_world2pix(world_coords, origin)
    print "  - Round trip matches:", np.allclose(pix_coords, pix_coords2)

    return wcs


def construct_wcs_from_scratch():
    print "Constructing WCS from a dictionary header:"
    header = {'NAXIS': 2,
              'NAXIS1': 2048,
              'CTYPE1': 'RA---TAN',
              'CRVAL1': 22.828128476,
              'CRPIX1': 1025.0,
              'CDELT1': 1.0,
              'NAXIS2': 1489,
              'CTYPE2': 'DEC--TAN',
              'CRVAL2': -0.945969070278,
              'CRPIX2': 745.0,
              'CDELT2': 1.0,
              'CD1_1': -1.06502217270E-08,
              'CD1_2': 0.000109979647614,
              'CD2_1': 0.000109949614203,
              'CD2_2': 3.21789868737E-09}
    wcs = WCS(header)
    
    # Test a round-trip from pix to world and back
    origin = 0
    pix_coords = np.arange(1, 11).reshape((5, 2))
    world_coords = wcs.wcs_pix2world(pix_coords, origin)
    pix_coords2 = wcs.wcs_world2pix(world_coords, origin)
    print "  - Round trip matches:", np.allclose(pix_coords, pix_coords2)

    return wcs


def construct_HPX_wcs():
    print "Constructing an HPX healpix projection"
    header = {'NAXIS': 2,
              'CTYPE1': 'RA---HPX',
              'CTYPE2': 'DEC--HPX',
             }
    wcs = WCS(header)
    
    # Test a round-trip from pix to world and back
    origin = 0
    pix_coords = np.arange(1, 11).reshape((5, 2))
    world_coords = wcs.wcs_pix2world(pix_coords, origin)
    pix_coords2 = wcs.wcs_world2pix(world_coords, origin)
    print "  - Round trip matches:", np.allclose(pix_coords, pix_coords2)

    return wcs


def plot_HPX_RA_DEC_GRID():
    header = {'NAXIS': 2,
              'CTYPE1': 'RA---HPX',
              'CTYPE2': 'DEC--HPX',
             }
    wcs = WCS(header)
    RA = np.linspace(0, 360, 40)
    DEC = np.linspace(-90, 90, 20)

    origin = 0
    RA, DEC = np.meshgrid(RA, DEC)

    RADEC = np.vstack([RA.ravel(), DEC.ravel()]).T
    
    x, y = wcs.wcs_world2pix(RADEC, origin).T


    plt.plot(x, y, '.k')
    plt.gca().set_aspect('equal')
    plt.show()


def plot_image_pixels_in_HPX():
    header = {'NAXIS': 2,
              'CTYPE1': 'RA---HPX',
              'CTYPE2': 'DEC--HPX',
             }
    hpx = WCS(header)

    hdulist = get_stripe82_file(301, 2700, 1, 'g', 18)
    img = WCS(hdulist[0].header)

    origin=0
    xpix, ypix = np.meshgrid(*[np.arange(0, hdulist[0].header[key], 30)
                               for key in ['NAXIS1', 'NAXIS2']])
    Xpix = np.vstack(map(np.ravel, (xpix, ypix))).T

    Xworld = img.wcs_pix2world(Xpix, origin)
    Xhpx = hpx.wcs_world2pix(Xworld, origin)

    x, y = Xhpx.T

    plt.plot(x, y, '.k')
    plt.gca().set_aspect('equal')
    plt.xlim(0, 360)
    plt.ylim(-90, 90)
    plt.show()


def plot_LSST_wcs(filename):
    hdulist = fits.open(filename)
    wcs = WCS(hdulist[1].header)

    hpx_header = {'NAXIS': 2,
                  'CTYPE1': 'RA---HPX',
                  'CTYPE2': 'DEC--HPX',
              }
    hpx = WCS(hpx_header)
    
    origin = 0
    xpix, ypix = np.meshgrid(*[np.arange(0, hdulist[1].header[key], 200)
                               for key in ['NAXIS1', 'NAXIS2']])
    Xpix = np.vstack(map(np.ravel, (xpix, ypix))).T

    for func in [wcs.wcs_pix2world, wcs.all_pix2world]:
        Xworld = func(Xpix, origin)
        Xhpx = hpx.wcs_world2pix(Xworld, origin)

        x, y = Xhpx.T

        plt.plot(x, y, '.')
        plt.gca().set_aspect('equal')
        plt.xlim(0, 360)
        plt.ylim(-90, 90)
    plt.show()


def test_SIP_distortion(filename):
    hdulist = fits.open(filename)
    wcs = WCS(hdulist[1].header)

    Xpix = [[1, 1]]
    Xfoc = wcs.sip_pix2foc(Xpix, 0)
    print wcs.sip_foc2pix(Xfoc, 0)
    return

    Xpix = 100 * np.random.random((40, 2))

    print "Conversion without SIP"
    Xworld = wcs.wcs_pix2world(Xpix, 0)
    Xpix2 = wcs.wcs_world2pix(Xworld, 0)
    print " - match:", np.allclose(Xpix, Xpix2)

    print "SIP by hand"
    Xworld = wcs.all_pix2world(Xpix, 0)
    tmp = wcs.sip_pix2foc(Xpix, 0)
    Xworld2 = wcs.wcs_pix2world(tmp, 0)
    print " - match:", np.allclose(Xworld, Xworld2)

    print "Conversion with SIP"
    tmp1 = wcs.sip_pix2foc(Xpix, 0)
    Xworld = wcs.wcs_pix2world(tmp1, 0)
    tmp2 = wcs.wcs_world2pix(Xworld, 0)
    Xpix2 = wcs.sip_foc2pix(tmp2, 0)
    print " - match:", np.allclose(tmp1, tmp2)
    print " - match:", np.allclose(Xpix, Xpix2)
    print Xpix - Xpix2


def extract_header(filename):
    header = fits.open(filename)[1].header
    import IPython; IPython.embed()
    
    
    
if __name__ == '__main__':
    filename = "/Users/jakevdp/research/LSST_IMGS/v865833781-fr/R21/S12.fits"
    extract_header(filename)
    #test_SIP_distortion(filename)
    
