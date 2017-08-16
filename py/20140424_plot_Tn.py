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
# Same setup as before, no coupler, +2.9 V on each amp
#d = bmxdata.BMXFile(datadir+'170424_1655.data')
#t290 = [1,90]
#t77 = [180,210]

# Now +3.0 V on each amp
#d = bmxdata.BMXFile(datadir+'170424_1712.data')
#t290 = [2,55]
#t77 = [140,160]

# Still +3 V, new LNA and 2nd stage amp
#d = bmxdata.BMXFile(datadir+'170424_1722.data')
#t290 = [45,50]
#t77 = [122,135]

# Now decrease ADC gain from 200 mV full ragen to 500 mV full range
#d = bmxdata.BMXFile(datadir+'170424_1730.data')
#t290 = [1,2]
#t77 = [3,4]

# Holy moly, you can't see the spectrum change at all! Now turn the amps off
# (V=0) and look at spectrum, it's basically zero. The time series goes between
# -1 and +1 ADU in places, but much less than with the amps on. Turn the amps
# on, see the spectrum go up and more hash in the time series. But amazing that
# you can't see the signal change. Points to either very high second stage amp
# noise or digitizer noise. 
#d = bmxdata.BMXFile(datadir+'170424_1737.data')
#t290 = [1,2]
#t77 = [3,4]

# With amps off, noise floor in power units is 1e11 at 500 mV full range, 1e11
# at 200 mV full range. Seems that the decrease is not enough, perhaps reaching
# some sort of ADC additive noise floor that does not care about the gain? Clear
# we probably need more gain.

# Add a second -3dB + LNA after second stage amp
d = bmxdata.BMXFile(datadir+'170424_1800.data')
t290 = [107,122]
t77 = [200,220]




# Get regions of constant temperature
t = arange(p.shape[0])*d.deltaT
ind290 = t2ind(t,t290)
ind77 = t2ind(t,t77)



p = d.data['chan1_0']
f = d.freq[0] + 1100


# Nominal freq index
ind = 800
f0 = f[ind]


# Plot

clf();
subplot(2,2,1)
plot(t,p[:,ind])
xlabel('time (s)')
ylabel(r'ADU$^2$ @ {:0.1f} MHz'.format(f0))
xlim(0,400)
for k in range(2):
    plot([t[ind290[k]],t[ind290[k]]],[0,2.5e14],'--k')
    plot([t[ind77[k]],t[ind77[k]]],[0,2.5e14],'--k')
    ylim(0,2.5e14)
    grid('on')

# Get average spectra
p290 = p[range(ind290[0],ind290[1])].mean(axis=0)
p77 =  p[range(ind77[0], ind77[1])].mean(axis=0)


subplot(2,2,2)
plot(f,p290,label='290 K avg.')
plot(f,p77,label='77 K avg.')
xlabel('f (MHz)')
ylabel(r'ADU$^2$')
legend()
ylim(0,2.5e14)
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
ylim(0,200)
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
