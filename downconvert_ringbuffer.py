#!/usr/bin/env python

import bmxdata
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft, irfft


#d = bmxdata.BMXRingbuffer('/astro/u/anze/bmxdata/raw/2011/ring/201107_225100_D1.ring')
d = bmxdata.BMXRingbuffer('/data/201107_225100_D1.ring')
del d.datad1c1
del d.datad1c2
del d.datad0c2

N = d.datad0c1.shape[0]
outdata = []
Nch=16
for chunk in range(1,Nch-1):
    i=(chunk-1)*N//Nch
    j=(chunk+2)*N//Nch
    print(i,j)
    T = (j-i)/1.1e9
    DeltaF = 1/T
    u = int(5e6/DeltaF)
    v = int(6e6/DeltaF)
    while ((v-u)%3!=1):
        v+=1
    data = d.datad0c1[i:j].astype(np.float)
    fft = rfft(data)
    print(len(data),len(fft),'xx',len(fft)*DeltaF)
    fft = fft[u:v]
    fft[0]=0
    fft[-1]=0.0
    data = irfft(fft)
    assert(len(data)%3==0)
    ip=len(data)//3
    jp=len(data)*2//3
    print (N/Nch, jp-ip)
    if chunk==2:
        np.savez('d1.npz',(data[jp-1024:jp+1024]))
    elif chunk==3:
        np.savez('d2.npz',(data[ip-1024:ip+1024]))
    print(data.shape)
    data = data[ip:jp]
    print(data.shape)
    outdata.append(np.copy(data))
    print (len(outdata), np.hstack(outdata).shape,'yy')

outdata = np.hstack(outdata)
print(outdata.shape)
np.savez('dconv.npz',outdata)
