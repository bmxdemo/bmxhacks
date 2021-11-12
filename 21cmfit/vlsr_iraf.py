import numpy as np
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy import coordinates as coord
import astropy.units as u
from astropy.time import Time

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

def vlsr(ra, dec, mjd):
    """RA, Dec in degrees. Return corrections that must be added to observed
    radial velocity to correct to local standard of rest.
    from http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?rvcorrect
    """

    dtor = np.pi / 180.
    rtod = 180 / np.pi

    # v_annual
    t = mjd / 36525.

    manom = 358.47583+t*(35999.04975-t*(0.000150+t*0.000003))
    oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503))
    lperi = 101.22083+t*(1.7191733+t*(0.000453+t*0.000003))
    eccen = 0.01675104-t*(0.00004180+t*0.000000126)
    
    tanom = manom + (2 * eccen - 0.25 * eccen**3) * np.sin(manom*dtor) \
            + 1.25 *eccen**2 * np.sin(2*manom*dtor) + \
            13./12. * eccen**3 * np.sin(3*manom*dtor)

    v = ((2*np.pi * 149598500.) / (365.2564 * 86400.)) / \
        np.sqrt(1. - eccen**2)
    

    slong = lperi + tanom + 180

    ep = 23.4 # obliquity of ecliptic
    b = rtod * np.arcsin(np.sin(dec*dtor)*np.cos(ep*dtor) - np.cos(dec*dtor)*np.sin(ra*dtor)*np.sin(ep*dtor))
    l = rtod * np.arccos(np.cos(dec*dtor)*np.cos(ra*dtor)/np.cos(b*dtor))

    v_annual = v * np.cos(b*dtor) * (np.sin((slong-l)*dtor) - eccen*np.sin((lperi-l)*dtor))


    # v_solar
    ra_apex = 18.0 * 360./24.
    dec_apex = 30.0
    
    v_sun = 19.5*(np.cos(ra_apex*dtor)*np.cos(dec_apex*dtor)*np.cos(ra*dtor)*np.cos(dec*dtor) + 
                  np.sin(ra_apex*dtor)*np.cos(dec_apex*dtor)*np.sin(ra*dtor)*np.cos(dec*dtor) +
                  np.sin(dec_apex*dtor)*np.sin(dec*dtor))

    # v_rot
    v_rot = 0.0



    return v_sun, v_annual, v_rot



    
