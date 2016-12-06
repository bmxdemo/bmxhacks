#!/usr/bin/env python
import bmxdata as bmx
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import pandas as pd

filelist="""161130_1657.tone.data
161130_1800.tone.data
161130_1900.tone.data
161130_2000.tone.data
161130_2100.tone.data
161130_2200.tone.data
161130_2300.tone.data
161201_0000.tone.data
161201_0100.tone.data
161201_0200.tone.data
161201_0300.tone.data
161201_0400.tone.data
161201_0500.tone.data
161201_0600.tone.data
161201_0700.tone.data
161201_0800.tone.data
161201_0900.tone.data
161201_1000.tone.data
161201_1100.tone.data
161201_1200.tone.data
161201_1300.tone.data
161201_1400.tone.data
161201_1500.tone.data
161201_1600.tone.data
161201_1700.tone.data
161201_1800.tone.data
161201_1900.tone.data
161201_2000.tone.data
161201_2100.tone.data
161201_2200.tone.data
161201_2300.tone.data
161202_0000.tone.data
161202_0100.tone.data
161202_0200.tone.data
161202_0300.tone.data
161202_0400.tone.data
161202_0500.tone.data
161202_0600.tone.data
161202_0700.tone.data
161202_0800.tone.data
161202_0900.tone.data
161202_1000.tone.data
161202_1100.tone.data
161202_1200.tone.data
161202_1300.tone.data
161202_1400.tone.data
161202_1500.tone.data
161202_1600.tone.data
161202_1700.tone.data
161202_1800.tone.data
161202_1900.tone.data
161202_2000.tone.data
161202_2100.tone.data
161202_2200.tone.data
161202_2300.tone.data
161203_0000.tone.data
161203_0100.tone.data
161203_0200.tone.data
161203_0300.tone.data
161203_0400.tone.data
161203_0500.tone.data
161203_0600.tone.data
161203_0700.tone.data
161203_0800.tone.data
161203_0900.tone.data
161203_1000.tone.data
161203_1100.tone.data
161203_1200.tone.data
161203_1300.tone.data
161203_1400.tone.data
161203_1500.tone.data
161203_1600.tone.data
161203_1700.tone.data
161203_1800.tone.data
161203_1900.tone.data
161203_2000.tone.data
161203_2100.tone.data
161203_2200.tone.data
161203_2300.tone.data
161204_0000.tone.data
161204_0100.tone.data
161204_0200.tone.data
161204_0300.tone.data
161204_0400.tone.data
161204_0500.tone.data
161204_0600.tone.data
161204_0700.tone.data
161204_0731.tone.data
161204_0900.tone.data
161204_1000.tone.data
161204_1100.tone.data
161204_1200.tone.data
161204_1300.tone.data
161204_1400.tone.data
161204_1500.tone.data
161204_1600.tone.data
161204_1700.tone.data
161204_1800.tone.data
161204_1900.tone.data
161204_2000.tone.data
161204_2100.tone.data
161204_2200.tone.data
161204_2300.tone.data
161205_0000.tone.data
161205_0100.tone.data
161205_0200.tone.data
161205_0300.tone.data
161205_0400.tone.data
161205_0500.tone.data
161205_0600.tone.data
161205_0700.tone.data
161205_0800.tone.data
161205_0900.tone.data
161205_1000.tone.data""".split("\n")


sw=0.0
#filelist=filelist[0:2]
for fname in filelist:
    d=bmx.BMXFile('../bmxdaq/data/'+fname)
    da=d.data['chan1_0']
    da=np.log(da)-35.14
    da[:,1638]=0.0
    w=len(da)
    cov=np.cov(da,rowvar=False,bias=True)
    cmean=da.mean(axis=0)
    cov=cov+np.outer(cmean,cmean)
    ###print (cov[0,0],cmean[0])
    ###print (cmean[:10])
    if sw==0:
        mean=cmean*w
        tcov=cov*w
        sw=w
    else:
        mean+=cmean*w        
        tcov+=cov*w
        sw+=w

mean/=sw
np.save('mean.dat',mean+35.14)
tcov/=sw
print (tcov[0,0])
tcov-=np.outer(mean,mean)
print (tcov[0,0],mean[0])

print(tcov[:4,:4])

print("doing eigen")
evl,evc=la.eig(tcov)
print(evl[:10])
np.save('evl.dat',evl)
np.save('evc.dat',evc)
print("done.")

"""

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



"""
