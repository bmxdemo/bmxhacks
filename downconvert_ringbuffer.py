#!/usr/bin/env python

import bmxdata
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, irfft
import glob
import os


def downbase(data, split=1024, N=2**11, pad=1):
    Nb=N*split
    Nch=len(data)//Nb
    print ("Nchunk=",Nch)
    fftsize = N*(2*pad+1)//2+1
    hambig = np.hamming((2+pad)*Nb)#*2+2)[1:-2:2]
    hamsmall = np.hamming((2+pad)*N)*split
    output = np.zeros((split,Nch-2*pad,N),np.float32)

    for i in range(pad,Nch-pad):
        print (i,'/',Nch-pad)
        subdata=data[(i-pad)*Nb:(i+2*pad)*Nb]*hambig

        fft = np.fft.rfft(subdata)

        fftws = np.zeros(fftsize, np.complex)
        for j in range(split):
            ofs = 1+j*(fftsize-1)
            fftws[1:] = fft[ofs:ofs+fftsize-1]
            chunk = np.fft.irfft(fftws)/hamsmall
            output[j,i-pad,:] = chunk[N:(1+pad)*N]
    return output




def process_ring(fn):
    outroot = "/astro/u/anze/bmxdata/baseband_ring/"+fn[fn.rfind('/')+1:fn.rfind('_D')]
    ndx = 4 if "D2" in fn else 0
    print (outroot)
    if not os.path.exists(outroot):
        os.mkdir(outroot)
    if os. path.exists(outroot+f"/{ndx}_0001.npy"):
        print(fn,"already processed. Skipping")
        return

    d = bmxdata.BMXRingbuffer(fn)
    
    for cc,dx in enumerate([d.datad0c1, d.datad0c2, d.datad1c1, d.datad1c2]):
        out = downbase(dx)
        nchan,ntime,chunksize = out.shape
        print (nchan,ntime,chunksize)
        for i in range(nchan):
            np.save(outroot+f'/{ndx+cc}_{i:04d}.npy',out[i,:,:])
       



for fn in glob.glob('/astro/u/anze/bmxdata/raw/2011/ring/*.ring'):
    process_ring(fn)

