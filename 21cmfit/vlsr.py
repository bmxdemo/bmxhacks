import numpy as np
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy import coordinates as coord
import astropy.units as u
from astropy.time import Time

telescope_loc = EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)

## Nobody knows where this was stole from
## But it was buggy -- I had to fix arccos for lam and lam_s

def vlsr(ra_deg, dec_deg, mjd):
    """RA, Dec in degrees. Return corrections that must be added to observed
    radial velocity to correct to local standard of rest."""

    ra_apex = (18.0 * 360./24.)*np.pi/180
    dec_apex = 30.0*np.pi/180
    
    ra = ra_deg * np.pi/180
    dec = dec_deg * np.pi/180

    # Solar motion w.r.t. LSR in km/s
    v_sun = 19.5*(np.cos(ra_apex)*np.cos(dec_apex)*np.cos(ra)*np.cos(dec) + 
                  np.sin(ra_apex)*np.cos(dec_apex)*np.sin(ra)*np.cos(dec) +
                  np.sin(dec_apex)*np.sin(dec))



    # Earth orbital motion w.r.t. sun
    c = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    ep = 23.4 * np.pi/180 # obliquity of ecliptic

    beta = np.arcsin(np.sin(dec)*np.cos(ep) - np.cos(dec)*np.sin(ra)*np.sin(ep))
    coslam = np.cos(dec)*np.cos(ra)/np.cos(beta)
    sinlam = np.cos(dec)*np.sin(ra)*np.cos(ep)+ np.sin(dec)*np.sin(ep)
    lam=np.arccos(coslam)
    lam[sinlam<0]*=-1
    

    time = Time(mjd, format='mjd')
    csun = coord.get_sun(time)
    ras  = csun.ra.deg  * np.pi/180
    decs = csun.dec.deg * np.pi/180
    beta_s = np.arcsin(np.sin(decs)*np.cos(ep) - np.cos(decs)*np.sin(ras)*np.sin(ep))
    coslam_s = np.cos(decs)*np.cos(ras)/np.cos(beta_s)
    sinlam_s = np.cos(decs)*np.sin(ras)*np.cos(ep)+ np.sin(decs)*np.sin(ep)
    lam_s = np.arccos(coslam_s)
    lam_s[sinlam_s<0]*=-1



    V = 29.974 # km/s
    e = 0.0167
    mjd0 = Time('1900-01-01', format='iso').mjd
    T = (mjd - mjd0)/365.25/100.
    gamma = 281.0 + 13./60 + 15.00/3600 + 6189.03*T/3600 + 1.63*T**2/3600
    gamma = gamma * np.pi/180

    v_earth = V*np.cos(beta)*(e*np.sin(gamma + lam) -
                 np.sin(lam - lam_s  ))


    # Earth rotational axis
    Ve = 0.465 # km/s
    h = 0 # assumes zenith pointing, hour angle always zero
    v_rot = Ve * np.sin(h)*np.cos(dec)*np.cos(telescope_loc.latitude.rad)
    
    return v_sun, v_earth, v_rot



    
