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
    t=0
    tsince=[]
    for i,r in enumerate(b.radarOn):
        if r:
            t=0
        else:
            t+=1
        tsince.append(t)
    tsince=np.array(tsince)
    chan='chan1_0'
    Ron=b.data[b.radarOn]['chan1_0'].mean(axis=0)
    Roff=[]

    st=3
    nst=98/st+1
    for i in range(nst):
        imin=i*st
        imax=(i+1)*st
        w=(tsince>=imin) & (tsince<imax)
        print i,w.sum()
        Roff.append(b.data[w]['chan1_0'].mean(axis=0))

    Rmin=np.array(Roff).min(axis=0)
    for sens in range(1):
        plt.figure()
        plt.plot(b.freq[0],Ron/Rmin-1,'g-',label='Radar on')
        for i in range(nst):
            plt.plot(b.freq[0],Roff[i]/Rmin-1,label='Dt=%3.2f s'%(0.122*(st*i+5)),color=(1.0*i/nst,0,1.-i*1.0/nst))
        plt.xlabel('MHz')
        plt.ylabel('Delta P')
        plt.legend(fontsize=10)
        #plt.ylim(0,2e16/10**sens)
        plt.tight_layout()
        plt.semilogy()
        plt.show()
        #plt.savefig("dRadar_%i.png"%sens)
        
    


if __name__=="__main__":
    main()
    
