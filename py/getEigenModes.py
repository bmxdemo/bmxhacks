#!/usr/bin/env python
import bmxdata as bmx
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import pandas as pd

d=bmx.BMXFile('data/161130_1800.tone.data')

mxf=[]
mx=[]
for line in d.data['chan1_1']:
    i=abs(line).argmax()
    mxf.append(d.freq[1][i])
    mx.append(line[max(0,i-20):i+20].sum()/3e18)
mx=np.array(mx)

mx=(mx/mx.mean()-1.0)

if False:
    d.nullBin(1638)
    da=d.data['chan1_0']
    da/=1e15

else:
    da=d.data['chan1_0']
    da=np.log(da)
    da[:,1638]=0.0

    
    
fr=d.freq[0]
N=len(da[0,:])
Nr=len(da)
print "N=",N
mean=da.mean(axis=0)
if False:
    cov=np.cov(da,rowvar=False)
    plt.figure()
    plt.plot(fr,mean)
    plt.show()
    plt.figure()
    plt.imshow(np.log(cov))
    plt.colorbar()
    plt.show()

if False:
    evl,evc=la.eig(cov)
    plt.figure()
    plt.hist(evl)
    plt.show()
    print evl[:10]
    plt.figure()
    np.save('evl.dat',evl)
    np.save('evc.dat',evc)
    for i in range(6):
        plt.subplot(3,2,i+1)
        plt.plot(fr,evc[:,i],label=str(i))
        plt.legend()
    plt.show()
evl=np.load('evl.dat.npy')
evc=np.load('evc.dat.npy')

print evl[0:5]
print evl.sum()
stop()


## now project eigen modes
da-=np.outer(np.ones(Nr),mean)

comp=[np.zeros(Nr) for i in range(5)]
for cc in range(5):
    Nc=24
    step=Nr/Nc
    for i in range (Nc):
        st=i*step
        en=(i+1)*step
        print i,st,en
        comp[cc][st:en]=np.dot(da[st:en,:],evc[:,cc])
if False:
    plt.figure()
    for i,l in enumerate(comp):
        plt.plot(l+0.1*i,label='component'+str(i))
        plt.xlabel('time stamp')
        plt.ylabel('component ampl')
    plt.legend()
    plt.show()


mx=pd.rolling_mean(mx,400)
c0=pd.rolling_mean(comp[0],400)
plt.figure()
#plt.plot(mx*(-40)-0.01,'r-', label='tone')
#plt.plot(c0,'b-', label='comp0')
plt.plot(mx,c0,'bo')
#plt.xlabel('time')
#plt.ylabel('smoothed on 40s running avg')
plt.legend()
plt.show()



