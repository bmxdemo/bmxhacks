#
# HI4PI for bmx calculator
#

import numpy as np
import fitsio
from astropy.wcs import WCS


class BMX4PI:
    raw_hi4pi='/astro/u/anze/bmxdata/HI4PI/CAR.fits'
    bmx_hi4pi='/astro/u/anze/bmxdata/HI4PI/bmxslice.fits'

    @staticmethod
    def create_slice (minvel=-400e3, maxvel=+500e3, mindec=28, maxdec=52):
        print ("Seetting up WCS.")
        w = WCS(BMX4PI.raw_hi4pi)
        print ("Reading %s ... "%BMX4PI.raw_hi4pi)
        da=fitsio.read(BMX4PI.raw_hi4pi)
        print ("Slicing...")
        NFR,NDEC,NRA=da.shape
        ve=w.wcs_pix2world(NRA/2,NDEC/2,np.arange(NFR),0)[2]
        ves=np.argmax(ve>=minvel)
        vee=np.argmax(ve>maxvel)
        de=w.wcs_pix2world(NRA/2,np.arange(NDEC),NFR/2,0)[1]
        des=np.argmax(de>=+mindec)
        dee=np.argmax(de>+maxdec)
        ra=w.wcs_pix2world(np.arange(NRA),NDEC/2,NFR/2,0)[0]
        rast=np.argmax(~np.isnan(ra))
        rae=len(ra)-np.argmax(~np.isnan(ra[::-1]))
        ra[ra>180.]-=360
        ra+=180
        vels=ve[ves:vee]
        decs=de[des:dee]
        ras=ra[rast:rae]
        skyslice=da[ves:vee,des:dee,rast:rae]
        print ("Dumping to %s ..."%BMX4PI.bmx_hi4pi)
        fitsio.write(BMX4PI.bmx_hi4pi,ras,clobber=True)
        fitsio.write(BMX4PI.bmx_hi4pi,decs)
        fitsio.write(BMX4PI.bmx_hi4pi,vels)
        fitsio.write(BMX4PI.bmx_hi4pi,skyslice)
        return (ras,decs,vels,skyslice)


    def __init__ (self, prepackaged=None):
        if prepackaged is not None:
            ras,decs,vels,skyslice=prepackaged
        else:
            print ("Reading %s ..."%BMX4PI.bmx_hi4pi)
            ra=fitsio.read(BMX4PI.bmx_hi4pi,ext=0)
            dec=fitsio.read(BMX4PI.bmx_hi4pi,ext=1)
            vel=fitsio.read(BMX4PI.bmx_hi4pi,ext=2)
            skyslice=fitsio.read(BMX4PI.bmx_hi4pi,ext=3)
        self.ra=ra
        self.dec=dec
        self.vel=vel
        self.sky=skyslice

    def Observed (self,ra,dec,sigmalat,sigmalon):
        """ Returns an observed signal for a given ra,dec fit
            with a certain sigma in beam in lat and longitude """
        pass
            
            
        
