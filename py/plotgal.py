#!/usr/bin/env python

import numpy as np
import sys
import bmxdata as bmx
import glob
import astropy.units as u
from scipy.interpolate import interp1d
from astropy.coordinates import EarthLocation, AltAz,SkyCoord
from astropy.time import Time
f21=1420.4

telescope_loc=EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m)
prepareData=True
day=31


def getB(da):

    mjdin=np.linspace(da.data['mjd'][0],da.data['mjd'][-1],100) ## good enough for our resolution
    time=Time(mjdin,format='mjd')
    point=AltAz(alt=90*u.deg, az=0*u.deg,location=telescope_loc,obstime=time)
    skyi=point.transform_to(SkyCoord(0*u.deg,0*u.deg,frame='galactic'))
    b=interp1d(mjdin,skyi.b.degree)(da.data['mjd'])
    return b

def main():
    f=sorted(glob.glob("/gpfs01/astro/workarea/bmxdata/1708%02d/1708%02d_1*.data"%(day-1,day-1))+
        glob.glob("/gpfs01/astro/workarea/bmxdata/1708%02d/1708%02d_2*.data"%(day-1,day-1))+
        glob.glob("/gpfs01/astro/workarea/bmxdata/1708%02d/1708%02d_0[0,1,2,3,4,5]*.data"%(day,day)))
    imgX=None
    bmin=-22
    bmax=30
    nB=(bmax-bmin)+1
    if prepareData:
        for fname in f:
            print fname
            bf=bmx.BMXFile(fname)
            bf.getRadar()
            galb=getB(bf)

            if imgX is None:
                imin,imax,fmin,fmax=bf.f2iminmax(1415.,1424.)
                nF=imax-imin
                imgX=np.zeros((nF,nB))
                imgY=np.zeros((nF,nB))
                sc=np.zeros(nB)
            for ib,b in enumerate(range(bmin,bmax+1)):
                #wh=np.where((not bf.radarOn) & (galb>b) & (galb<b+1))
                wh=np.where((galb>b) & (galb<b+1))
                if len(wh[0])>0:
                    print b
                    imgX[:,ib]+=(bf.data[wh]['chan1_0'])[:,imin:imax].sum(axis=0)
                    imgY[:,ib]+=(bf.data[wh]['chan2_0'])[:,imin:imax].sum(axis=0)
                    sc[ib]+=len(wh[0])
        imgX/=np.outer(np.ones(nF),sc+1e-3)
        imgY/=np.outer(np.ones(nF),sc+1e-3)
        np.save('g'+str(day),(imgX,imgY))

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    imgX,imgY=np.load("g%i.npy"%day)
    for i,img in enumerate([imgX,imgY]):
        ## fft filter
        print "FFT"
        if False:
            t=np.fft.rfft2(img)
            t[:,30:]=0 # fast varying in b
            t[:,:2]=0 # non-varying in b

            t[:4,:]=0 #slow vary in f
            t[-4:,:]=0 #slowvary in f

            print t.shape
            img=np.fft.irfft2(t)
        
        plt.subplot(2,1,i+1)
        plt.imshow(img.T,interpolation="nearest" , aspect = "auto", extent=[1415-f21, 1424-f21, bmin,bmax],origin='lower')


        plt.xlabel("$\Delta$ MHz")
        plt.ylabel("galactic b")
        plt.colorbar()
    plt.show()



if __name__=="__main__":
    main()
    
        
    
