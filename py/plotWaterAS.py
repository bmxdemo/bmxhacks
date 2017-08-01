#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import bmxdata as bmx


def main():
    f = []
    f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170718_2000.data")
    if True:
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170718_2100.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170718_2200.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170718_2300.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0000.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0100.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0200.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0300.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0400.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0500.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0600.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0700.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0800.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_0900.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_1000.data")
        f.append("/gpfs01/astro/workarea/anze/bmxdata/170719_first/flight_rfi_20_4_170719_1100.data")
    bigar=[None,None]
    #fmin,fmax=1400., 1440.
    fmin,fmax=None,None
    #fnormmin,fnormmax=1320.,1340.
    fnormmin,fnormmax=None,None
    binSize=4
    cut=0
    tavg=600
    minmax=0.3
    for fname in f:
        b=bmx.BMXFile(fname)
        if False:
            b.plotWaterfall(nsamples=300,binSize=4)
            stop()
        b.getRadar()
        b.filterRadar()
        if fnormmin is not None:
            b.normalizeOnRegion(fnormmin,fnormmax)
        imin,imax,fminp,fmaxp=b.f2iminmax(fmin,fmax,binSize=binSize)
        for n in range(2):
           arr = []
           for i in range(b.nSamples):
               arr.append(b.data[i]['chan' + str(n+1)+'_' + str(cut)][imin:imax])
               arr[i] = np.reshape(arr[i],(-1, binSize )) #bin frequencies into groups of 4
               arr[i] = np.mean(arr[i], axis = 1)  #average the 4
   	   arr=np.array(arr)
           if bigar[n] is None:
               bigar[n]=arr
           else:
               bigar[n]=np.vstack((bigar[n],arr))
    if tavg>1:
        nsize=len(bigar[0])/tavg
        for ch in [0,1]:
            for i in range(nsize):
                bigar[ch][i,:]=bigar[ch][tavg*i:tavg*(i+1),:].mean(axis=0)
            bigar[ch]=bigar[ch][:nsize,:]
    print '---------'        
    for ch in [0,1]:
        print bigar[ch].shape    
        mean=bigar[ch].mean(axis=0)
        if True:
            for i in range(len(bigar[ch])):
                bigar[ch][i]/=mean
            bigar[ch]-=1.
        ## fft filter
        print "FFT"
        if True:
            t=np.fft.rfft2(bigar[ch])
            t[5:-5,:]=0
            t[:,:50]=0
            bigar[ch]=np.fft.irfft2(t)
        
        
        plt.subplot(2,1,ch+1)
        if True:
            plt.imshow(bigar[ch],interpolation="nearest" , aspect = "auto", extent=[fminp, fmaxp, len(bigar[0])*.122016*tavg/3600, 0]
                       ,vmin=-minmax,vmax=+minmax
            )
        else:
            plt.imshow(bigar[ch],norm=colors.LogNorm(),interpolation="nearest" , aspect = "auto", extent=[fminp, fmaxp, len(bigar[0])*.122016*tavg/3600, 0]
                   )

        plt.xlabel("MHz")
        plt.ylabel("t[hours] since 8pm")
        plt.colorbar()
    plt.show()

            
            
        
        #b.plotWaterfall(nsamples=200,subtractMean=True,minmax=1)

        


if __name__=="__main__":
    main()
    
        
    
