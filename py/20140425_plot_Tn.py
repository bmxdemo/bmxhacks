import bmxdata
from numpy import *
from matplotlib.pyplot import *
ion()

def t2ind(t,t0):
    ind = []
    for k,val in enumerate(t0):
        ind.append(where(abs(t-val) == np.min(abs(t-val)))[0])
    return ind

# Get cooldown data
datadir = '/direct/astro+u/csheehy/python/bmxdaq/data/'


#############
# LNA1 SN: F371501603
# LNA2 SN: F556901639
# AMP2 SN: F566301633

# From last time...
# Same setup as before, no coupler, +2.9 V on each amp
# LNA -> -3dB -> filter -> AMP2 -> lowpass -> -30 dB
#d = bmxdata.BMXFile(datadir+'170424_1655.data')
#t290 = [1,90]
#t77 = [180,210]

## Setup is now: (fig1)
## LNA? -> lowpass -> amp2 -> bandpass -> ADC
## LNA @ +3V, amp2 @ +2.8V
#d = bmxdata.BMXFile(datadir+'170425_1545.data')
#t290 = [100,150]
#t77 = [250,260]

## Setup is now: (fig2)
## LNA1 -> -3dB -> LNA2 -> lowpass -> amp2 -> bandpass -> ADC
## LNA @ +3V, amp2 @ +2.8V
#d = bmxdata.BMXFile(datadir+'170425_1606.data')
#t290 = [20,60]
#t77 = [140,155]

## Setup is now:
## LNA2 -> -3dB -> LNA1 -> lowpass -> amp2 -> bandpass -> ADC
## LNA @ +3V, amp2 @ +2.8V
#d = bmxdata.BMXFile(datadir+'170425_1629.data')
#t290 = [40,50]
#t77 = [140,150]

## Maintain setup but ensure short coax length worse end is facing down
## (e.g. terminated) 
#d = bmxdata.BMXFile(datadir+'170425_1647.data')
#t290 = [10,20]
#t77 = [90,110]

## Flip coax
#d = bmxdata.BMXFile(datadir+'170425_1657.data')
#t290 = [50,100]
#t77 = [180,200]

## Now LNA2 -> -3dB -> amp2 (bent pin) -> lowpass -> 4x -3dB -> amp2 -> bandpass
## -> ADC
#d = bmxdata.BMXFile(datadir+'170425_1927.data')
#t290 = [20,50]
#t77 = [140,160]

## Setup is now: (fig3)
## LNA2 -> LNA1 -> lowpass -> amp2 -> bandpass -> ADC
## LNA @ +3V, amp2 @ +2.8V
#d = bmxdata.BMXFile(datadir+'170425_1942.data')
#t290 = [20,34]
#t77 = [85,95]

## Same as above but substitute very short terminator link for coax
#d = bmxdata.BMXFile(datadir+'170425_1949.data')
#t290 = [10,40]
#t77 = [110,120]

## Back to
# LNA2 -> LNA1 -> lowpass -> amp2 -> bandpass -> ADC
# But now use same voltage on all amps, +2.9 V
#d = bmxdata.BMXFile(datadir+'170425_2006.data')
#t290 = [20,35]
#t77 = [90,100]

## Identical to short terminator link but with crappy paper vapor shield to try
#to keep LNA warm
#d = bmxdata.BMXFile(datadir+'170425_2018.data')
#t290 = [20,40]
#t77 = [100,110]

p = d.data['chan1_0']
f = d.freq[0] + 1100


# Get regions of constant temperature
t = arange(p.shape[0])*d.deltaT
ind290 = t2ind(t,t290)
ind77 = t2ind(t,t77)


# Nominal freq index
ind = 800
f0 = f[ind]


# Plot

clf()
subplot(2,2,1)
plot(t,p[:,ind])
xlabel('time (s)')
ylabel(r'ADU$^2$ @ {:0.1f} MHz'.format(f0))
for k in range(2):
    yl=ylim()
    plot([t[ind290[k]],t[ind290[k]]],[0,yl[1]],'--k')
    plot([t[ind77[k]],t[ind77[k]]],[0,yl[1]],'--k')
grid('on')

# Get average spectra
p290 = p[range(ind290[0],ind290[1])].mean(axis=0)
p77 =  p[range(ind77[0], ind77[1])].mean(axis=0)


subplot(2,2,2)
semilogy(f,p290,label='290 K avg.')
semilogy(f,p77,label='77 K avg.')
xlabel('f (MHz)')
ylabel(r'ADU$^2$')
legend(loc='lower left')
#ylim(0,2.5e14)
grid('on')


# Compute noise temperature
T1 = 77.0
T2 = 290.0
P1 = p77
P2 = p290

# P = gT + n
g = (P2 - P1)/(T2-T1) # ADU^2 / K
n = P1 - g*T1
Tn = n/g # K

subplot(2,2,3)
plot(f, Tn)
xlabel('f (MHz)')
ylabel(r'$T_n$ (K)')
ylim(0,50)
grid('on')


# Calibrate all the spectra
pcal = zeros(p.shape)
for k,val in enumerate(p):
    pcal[k,:] = (val-n)/g
subplot(2,2,4)
plot(f,pcal[10:3600:50].T)
ylim(0,300)
xlabel('f (MHz)')
ylabel('T (K)')
