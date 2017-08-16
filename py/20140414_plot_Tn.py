import bmxdata
from numpy import *
from matplotlib.pyplot import *
ion()

# Get cooldown data
for kk in range(2):
    if kk==0:
        d = bmxdata.BMXFile('data/170414_1618.data') # with coupler
        ind290 = [20,850]; ind77  = [1320, 1500]
    else:
        d = bmxdata.BMXFile('data/170414_1628.data') # no coupler
        ind290 = [20,400]; ind77  = [920, 1045]

    p = d.data['chan1_0']
    f = d.freq[0] + 1100

    # Time in seconds
    t = arange(p.shape[0])*d.deltaT

    # Nominal freq index
    ind = 800
    f0 = f[ind]


    # Plot


    if kk==0:
        clf();
        subplot(2,2,1)
        plot(t,p[:,ind])
        xlabel('time (s)')
        ylabel(r'ADU$^2$ @ {:0.1f} MHz'.format(f0))
        xlim(0,400)
        for k in range(2):
            plot([t[ind290[k]],t[ind290[k]]],[0,1.2e13],'--k')
            plot([t[ind77[k]],t[ind77[k]]],[0,1.2e13],'--k')
            ylim(0,1.2e13)
            grid('on')

    # Get average spectra
    p290 = p[range(ind290[0],ind290[1])].mean(axis=0)
    p77 =  p[range(ind77[0], ind77[1])].mean(axis=0)

    if kk==0:
        subplot(2,2,2)
        plot(f,p290,label='290 K avg.')
        plot(f,p77,label='77 K avg.')
        xlabel('f (MHz)')
        ylabel(r'ADU$^2$')
        legend()
        ylim(0,1.5e13)
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
    if kk==0:
        lb = 'with coupler'
    else:
        lb = 'no coupler'
    plot(f, Tn, label=lb)
    xlabel('f (MHz)')
    ylabel(r'$T_n$ (K)')
    ylim(0,200)
    if kk==1:
        legend()
    grid('on')


    # Calibrate all the spectra
    if kk==0:
        pcal = zeros(p.shape)
        for k,val in enumerate(p):
            pcal[k,:] = (val-n)/g
        subplot(2,2,4)
        plot(f,pcal[10:3600:50].T)
        ylim(0,300)
        xlabel('f (MHz)')
        ylabel('T (K)')
