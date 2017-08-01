#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import bmxdata as bmx


def main():
    fname="/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0300.data"
    b=bmx.BMXFile(fname)
    b.getRadar()
    mspec=b.data[np.logical_not(b.radarOn)]['chan1_0'].mean(axis=0)
    mspec2=b.data[np.logical_not(b.radarOn)]['chan2_0'].mean(axis=0)
    mspecR=b.data[np.logical_not(b.radarOn)]['chanXR_0'].mean(axis=0)
    mspecI=b.data[np.logical_not(b.radarOn)]['chanXI_0'].mean(axis=0)
    mspecRr=b.data[b.radarOn]['chanXR_0'].mean(axis=0)
    mspecIr=b.data[b.radarOn]['chanXI_0'].mean(axis=0)


    for i in range(4):
        plt.subplot(4,1,i+1)
        imin=1024*i
        imax=imin+1024
        plt.plot(b.freq[0][imin:imax],mspec[imin:imax])
        plt.plot(b.freq[0][imin:imax],mspec2[imin:imax])
        plt.xlabel("MHz")
        plt.ylabel("P")
        plt.semilogy()
    plt.show()

    for i in range(4):
        plt.subplot(4,1,i+1)
        imin=1024*i
        imax=imin+1024
        plt.plot(b.freq[0][imin:imax],mspecR[imin:imax])
        plt.plot(b.freq[0][imin:imax],mspecI[imin:imax])
        plt.xlabel("MHz")
        plt.ylabel("P")
    plt.show()



    
    return









    Nfft=4096*4
    da=np.zeros(Nfft,complex)
    da[4096*2:4096*3]=mspec
    #da[:4096]=1.0+0j
    xi=np.fft.irfft(da)[:Nfft]
    dt=3e8/3.3e9
    me=dt*np.arange(len(xi))
    plt.plot(me,xi)
    plt.semilogx()
    plt.show()
            
    


if __name__=="__main__":
    main()
    
