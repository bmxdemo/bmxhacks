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
datadir = './data/'

g = []
n = []

# Channel X, Y
for j in range(2):
    if j==0:
        figure(1)
        clf()
        'Xpol, channel 1'
        d = bmxdata.BMXFile(datadir+'170816_1502.data')
        p = d.data['chan1_0']
        t290 = [60,80]
        t77 = [300,310]
    elif j==1:
        figure(2)
        clf()
        'Ypol, channel 2'
        d = bmxdata.BMXFile(datadir+'170816_1538.data')
        p = d.data['chan2_0']
        t290 = [150,160]
        t77 = [350,360]


    f = d.freq[0] + 1100

    # Get regions of constant temperature
    t = arange(p.shape[0])*d.deltaT
    ind290 = t2ind(t,t290)
    ind77 = t2ind(t,t77)


    # Example freq index for cooldown time series, choosen basically at random. 
    ind = 800
    f0 = f[ind]


    # Plot

    subplot(2,2,1)
    plot(t,p[:,ind])
    xlabel('time (s)')
    ylabel(r'ADU$^2$ @ {:0.1f} MHz'.format(f0))
    xlim(0,400)
    for k in range(2):
        plot([t[ind290[k]],t[ind290[k]]],[0,7e14],'--k')
        plot([t[ind77[k]],t[ind77[k]]],[0,7e14],'--k')
        ylim(0,7e14)
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
    ylim(0,7e14)
    grid('on')


    # Compute noise temperature
    T1 = 77.0
    T2 = 303.0
    P1 = p77
    P2 = p290

    # P = gT + n
    g.append( (P2 - P1)/(T2-T1) )# ADU^2 / K
    n.append( (P1 - g[j]*T1) )
    Tn = n[j]/g[j] # K

    subplot(2,2,3)
    plot(f, Tn)
    xlabel('f (MHz)')
    ylabel(r'$T_n$ (K)')
    ylim(0,50)
    grid('on')


    # Calibrate all the spectra
    pcal = zeros(p.shape)
    for k,val in enumerate(p):
        pcal[k,:] = (val-n[j])/g[j]
    subplot(2,2,4)
    plot(f,pcal[10:3600:50].T)
    ylim(0,300)
    xlabel('f (MHz)')
    ylabel('T (K)')
