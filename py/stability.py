#!/usr/bin/env python
import bmxdata as bmx
import matplotlib.pyplot as plt
import numpy as np
import sys
from numpy.fft import rfft

a1=[]
a2=[]
b1=[]
b2=[]
for fn in sys.argv[1:]:
    print fn
    d=bmx.BMXFile(fn)
    a1.append(d.getToneAmplFreq('chan1_0')[1])
    b1.append(d.getToneAmplFreq('chan1_1')[1])
    a2.append(d.getToneAmplFreq('chan2_0')[1])
    b2.append(d.getToneAmplFreq('chan2_1')[1])

print [x.shape for x in a1]
a1=np.hstack(a1)
a2=np.hstack(a2)
b1=np.hstack(b1)
b2=np.hstack(b2)


def proc (a,b,avg=1):
    rat=a/b
    rat/=rat.mean()
    rat-=rat.mean()
    if (avg>1):
        m=np.array([rat[i*avg:(i+1)*avg].mean() for i in range(len(rat)/avg)])
        return m
    else:
        return rat

def procx (a,b,a2,b2,avg=1):
    rat=a/b/(a2/b2)
    rat/=rat.mean()
    rat-=rat.mean()
    if (avg>1):
        m=np.array([rat[i*avg:(i+1)*avg].mean() for i in range(len(rat)/avg)])
        return m
    else:
        return rat

    
rat1=proc(a1,b1)
rat2=proc(a2,b2)
rat=0.5*(rat1+rat2)
ratx=procx(a1,b1,a2,b2)
a1l=proc(a1,np.ones(len(a1)))
a12l=proc(a1,a2)

plt.plot(a1l)
a1l-=ratx * (a1l*ratx).sum()/(ratx**2).sum()
a12l-=ratx * (a12l*ratx).sum()/(ratx**2).sum()
rat1-=ratx * (rat1*ratx).sum()/(ratx**2).sum()




def plotps (ar,l):
    ft=rfft(ar)
    pk=abs(ft**2)
    freq=(np.arange(len(pk)))/(0.122*len(ar)/3600.)
    ## log binning
    st,step=0,1.
    lpk=[]
    lfr=[]
    fa=1.02
    while True:
        if (st>len(pk)):
            break
        lpk.append(pk[st:st+int(step)].mean())
        lfr.append(freq[st:st+int(step)].mean())
        st+=int(step)
        step*=fa
    plt.plot(lfr,lpk,lw=2, label=l)

def plotax(ar,l):
    avg=100
    m=np.array([ar[i*avg:(i+1)*avg].mean() for i in range(len(rat)/avg)])
    plt.plot(m,label=l)
    
plt.figure()
plotax(a1l,'T1')
plotax(a12l,'T1 cal off CH2')
plotax(rat1,'T1 cal off T2')
plotax(ratx,'T1 cal off T2,Ch2')
plt.xlabel("sample")
plt.ylabel("rat")
plt.legend()
plt.tight_layout()
plt.savefig("pltr.pdf")

plt.figure()
plotps(a1l,'T1')
plotps(a12l,'T1 cal of CH2')
plotps(rat,'T1 cal off T2')
plotps(ratx,'T1 cal off T2,Ch2')
plt.xlabel("f[1/hr]")
plt.ylabel("P(f)")
plt.ylim(1e-5,1e5)
plt.legend()
plt.tight_layout()


plt.loglog()
plt.savefig("pltps.pdf")
#plt.show()

#plt.plot(rat2)




