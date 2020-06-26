# imports
import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

BMXDATADIR='/gpfs01/astro/workarea/bmxdata/reduced/fra'


def droneDataPair(FLY, mjdTweak=0, verbose=False, debug=False):
        """ associate BMX autospectra with DJI coordinate data 
            inputs
            FLY: str, FLY number from Yale dataset
            mjdTweak: float, offset of mjdBMX 
            return:
            outdict: dict with keys 
                'FLY'   ID of DJI flight #
                'BMXdir'data directory in /gpfs01/astro/workarea/bmxdata/fra
                'pol'   polarization direction
                'attn'  attenuation
                 mjd    mjd according to BMX
                 x, y, z        drone coordinates (in m) interpolated to BMX
                 autos  list of 4 arrays with waterfall data for the dishes sensitive to the polarization direction
            
                Datasets for 200312 Yale data:

                BMX     FLY     POL     DIR     ATT     STRT    UTC
                1434    340                     30      10:32   14:32
                1702    342     N       E       30      1.04    17:04   # + NS pass same FLY?
                1813    343     N       N       30      1.38    17:38   # WRONG
                1732    343     N       N       30      1.38    17:38   # CORRECTED
                2118    352     E       E       30      5.18    21:18
                2152    353     E       N       30      5.52    21:52

                1916    344     N       E       43      2.12    18:12   # WRONG
                1813    344     N       E       43      2.12    18:12   # CORRECTED
                1944    347     N       N       43      3.45    19:45
                2011    349     E       E       43      4.10    20:10   # BMX mjd != BMX autosp.shape[0]
                2047    351     E       N       43      4.50    20:50   """

                # match BMX and DJI data sources (DJI urls point to csv files on gdrive)
        if FLY=='340':
                BMXdir = '200312_1434_yale'
                pol = 'N'
                attn = 30
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vTXNtVhOmLK3kBG2QGg8F7NRFulz_2tgEioanTionUeYFHoDfKaXLznOodTlWQEUoeeHYCOXeD7wvgw/pub?output=csv'

        if FLY=='342':
                BMXdir = '200312_1702_yale'
                pol = 'N'
                attn = 30
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vQgJV-BZoWBHL7BkNWOviln2ivpDD06Cf8i2H0wZY3ZzDXRcoI30bo-OHHRG_9HKnNcI1nJwNG01RWg/pub?output=csv'
        if FLY=='343':
                BMXdir = '200312_1732_yale'
                pol = 'N'
                attn = 30
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vQI17JixrkQhALHEhW4XDLl-Vj6WBzWtilP_AV-nnN9IuIdVIb5wSCMZpR6yh-FzXBFKuOGJeuAfhgY/pub?output=csv'

        if FLY=='352':
                BMXdir = '200312_2118_yale'
                pol = 'E'
                attn = 30
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vRGmjldQZwSKIf4lL03wBm6mhHEiWDwIvu9cLQ_S4xB-y7RrAqGUfWuAYkF1d1XOzBGjgKnHl51MUB8/pub?output=csv'

        if FLY=='353':
                BMXdir = '200312_2152_yale'
                pol = 'E'
                attn = 30
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vS9GK5cSNfT8ihD8I3_DosLsCkLa19CK6lb7Xwm4iqHqlL5G_dRmErS6WNEJcjUjYsMwHmLZyBhcXoA/pub?output=csv'

        if FLY=='344':
                BMXdir = '200312_1813_yale'
                pol = 'N'
                attn = 43
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vRu66U2XG2pU5DAkAN49ZhbhJtP_n4UqoVBO3toJ5OTPJixO1cmudMVCrkvIFIpIrJAKjGvJpz4hKcW/pub?output=csv'

        if FLY=='347':
                BMXdir = '200312_1944_yale'
                pol = 'N'
                attn = 43
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vStDvMGZE0GDh31U3h4oaTykl_rO9O6k0CnR2Cd7qrYoyyZOzL1umiTiApoTMbYnJFc3EGRYhgkG5gy/pub?output=csv'

        if FLY=='349':
                BMXdir = '200312_2011_yale'
                pol = 'E'
                attn = 43
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vRBwdk_z86hsQwMiWIQAc0rNceAUV-LhKZC0wRmP1-jvYp_h4idlL9bMSWlWJW_FgoY0FW7XsqDTz7z/pub?output=csv'

        if FLY=='351':
                BMXdir = '200312_2047_yale'
                pol = 'E'
                attn = 43
                DJIurl = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vTplEvCskI_GO3yb7YDNFm8thTfK4mmoYpSh-s8anNRqQ5Yjg2ZLCYvEKl22uwz1FHTd60M2Lx9vDrR/pub?output=csv'

        # get BMX mjd and optionally tweak:
        mjdB = fits.open(BMXDATADIR +os.sep + BMXdir + os.sep + 'mjd.fits')[1].data.field(0) + mjdTweak

        # get BMX autospectra for the channels responding to drone polarization
        autos = []
        for chan in range(4):
                if pol=='N':
                        fn = BMXDATADIR + os.sep + BMXdir + os.sep + 'cut0' + os.sep + 'auto_' + str(chan + 1) + '.fits'
                        autos.append(fits.open(fn)[0].data)
                if pol=='E':
                        fn = BMXDATADIR + os.sep + BMXdir + os.sep + 'cut0' + os.sep + 'auto_' + str(chan+5) + '.fits'
                        autos.append(fits.open(fn)[0].data)



        # read the DJI csv file and convert to dict
        if (verbose): 
                print "reading csv file..."
        DJIdict = pd.read_csv(DJIurl).to_dict(orient='list')


        # convert lat, lon in degrees to x,y in m
        latMid = median(DJIdict['Lat'])
        m_per_deg_lat = 111132.954 - 559.822*cos(radians(2*latMid)) + 1.175*cos(radians(4*latMid))
        m_per_deg_lon = 111132.954*cos(radians(latMid))
        xp = array(DJIdict['Lon'])*m_per_deg_lon
        yp = array(DJIdict['Lat'])*m_per_deg_lat
        zp = array(DJIdict['hmsl'])     # already in m
        if (debug):
                print "xp,yp,zp ranges: ", xp.max()-xp.min(), yp.max()-yp.min(), zp.max()-zp.min()

        # convert datetimestamp to mjd
        # note datetimestamp only updates once per second
        # so get MJD of first timestamp and increment by (sec_elapsed - sec_elapsed[0])/86400
        t = Time(DJIdict['datetimestamp'])
        MJD0 = t.mjd[0]
        se = array(DJIdict['sec_elapsed']) - DJIdict['sec_elapsed'][0]
        mjdD = MJD0 + se/86400.

        # interpolate DJI coordinates to mjdB times:
        y = np.interp(mjdB, mjdD, yp)
        x = np.interp(mjdB, mjdD, xp)
        z = np.interp(mjdB, mjdD, zp)
        if (debug):
                print "x,y,z ranges ", x.max()-x.min(), y.max()-y.min(), z.max()-z.min()

        # coordinates of dish centers (from FLY338)
        # NESW by y,x,z
        dishLoc = ([  4.53866934e+06,  -6.12354322e+06,   1.70000000e+01],
        [  4.53866568e+06,  -6.12353927e+06,   1.70000000e+01],
        [  4.53866241e+06,  -6.12354284e+06,   1.70000000e+01],
        [  4.53866555e+06,  -6.12354663e+06,   1.70000000e+01])

        # compute angles (x=lon, y=lat) from each dish
        thetaX = []
        thetaY = []
        for chan in range(4):
	        thetaX.append(np.arctan((x - dishLoc[chan][1])/(z - dishLoc[chan][2])))
	        thetaY.append(np.arctan((y - dishLoc[chan][0])/(z - dishLoc[chan][2])))
        thetaX = np.array(thetaX)
        thetaY = np.array(thetaY)

        if (verbose):
                print "creating output dictionary and terminating"

        return {"FLY" : FLY, "BMXdir" :  BMXdir, "pol" : pol, "attn" : attn,  "autos" : autos, "mjdD" : mjdD, "mjdB" : mjdB, "x" : x, "y" : y, "z" : z, "thetaX" : thetaX, "thetaY" : thetaY} 
