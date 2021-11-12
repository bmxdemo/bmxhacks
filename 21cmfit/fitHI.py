import sys
sys.path+=['.','..','bmxreduce']


from mapmanager import mapmanager
from numpy import *
import healpy as hp
import reduce_coadd
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy import optimize as opt
import cPickle as cP
from glob import glob
from sklearn.linear_model import Ridge
import vlsr
import vlsr_iraf
from xkeckhelio import x_keckhelio
from matplotlib.pyplot import *
import numpy as np
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", dest="i", type="int", default=0)
(o, args) = parser.parse_args()

loaddata = True

if loaddata:

    #dind = [0,8,9,10,15,16,17,18,19,20,21,22,30,35,36,37,52,91,92,93,94]
    dind = [-1]

    fn = sort(glob('maps/bmx/real/galcross/*.npz'))
    r=reduce_coadd.coaddbyday([fn[dind[o.i]]],dodeglitch=False)

    f = r.f

    df = f[1]-f[0]
    be_lo = f-df/2
    be_hi = f+df/2

loadtemplates = False
if loadtemplates:

    fhi = loadtxt('maps/HI4PI/stitched/nuTable.txt')
    fhi = fhi[200:700] # must match glob range below

    fid = open('personal/HItemplates/galcross_skycoords.pickle','rb')
    c = cP.load(fid)
    fid.close()

    fwhmarrtemp = np.arange(3.5, 5.5, 0.5)
    darrtemp = np.arange(-10.0, -8.25, 0.25)
    temp0 = np.zeros((len(fwhmarrtemp), len(darrtemp), len(fhi), len(r.ra)))

    for di,d in enumerate(darrtemp):
        print(d)

        for fi, fwhm in enumerate(fwhmarrtemp):

            fn = sort(glob('personal/HItemplates/template_[2-6]??_fwhm{:0.2f}_d{:0.1f}.npy'.format(fwhm,d)))
            for kk, val in enumerate(fn):
                y = np.load(val)
                temp0[fi, di, kk, :] = np.interp(r.ra,c.ra,y)

    temp0 = temp0 / 2. # divide by 2 for single polarization???!!! CHECK!!!

    np.savez('personal/HItempload.npz',fhi=fhi, fwhmarrtemp=fwhmarrtemp, darrtemp=darrtemp, temp0=temp0)
else:
    x=np.load('personal/HItempload.npz')
    fhi=x['fhi']
    fwhmarrtemp=x['fwhmarrtemp']
    darrtemp=x['darrtemp']
    temp0 = x['temp0']

barr = [969,970,971,972,973]

#fwhmarr = np.arange(3.5,5.5,0.5)
#darr = np.arange(-10.0,-8.25,0.25)
#dtarr = np.arange(-2,1)
#dfarr = np.arange(-0.03, 0, np.abs(fhi[2]-fhi[1]))

fwhmarr = [5.0]
darr = [-8.75]
dtarr = [0]
dfarr = [-0.01]

chi2 = np.zeros((len(fwhmarr),len(darr),len(dtarr),len(dfarr)))
mod = None
coef = None
rr=Ridge(alpha=0, normalize=False, fit_intercept=False)



for i,fwhm in enumerate(fwhmarr):

    print ('{:d} of {:d}'.format(i+1,len(fwhmarr)))
    fwhmind = np.where(fwhmarrtemp == fwhm)[0][0]

    for j, d in enumerate(darr):

        #print ('{:d} of {:d}'.format(j+1,len(darr)))
        dind = np.where(darrtemp == d)[0][0]

        for k, dt in enumerate(dtarr):

            for l, df in enumerate(dfarr):

                aa = None
                for m,b in enumerate(barr):

                    # Data to fit
                    v = (r.T+r.modcpm+r.mod)[0,:,b]
                    v0 = (r.T+r.modcpm+r.mod)[0,:,967]

                    # Subtract off zero point shift
                    v = v-v0
                    #v = v - v[-1]

                    # Velocity correction. This must be inside the loop! 
                    dra = (r.ra[2]-r.ra[1])*dt
                    ra0 = r.ra + dt
                    dec0 = r.dec + d

                    dv0 = np.squeeze(np.array([x_keckhelio(ra0[u],dec0[u],r.mjd[u]) for u in range(r.ra.size)]))
                    dv1 = vlsr.vlsr(ra0, dec0, r.mjd)
                    dv = dv0 - dv1[0]

                    #dv = np.zeros_like(ra0)

                    z = dv / 3.0e5

                    temp = np.zeros_like(r.mjd)
                    for kk,val in enumerate(z):
                        #fobs = fhi/(1+val)
                        fobs = fhi*(1-val)
                        doind = (fobs>(be_lo[b]+df)) & (fobs<(be_hi[b]+df))
                        temp[kk] = temp0[fwhmind, dind, doind, kk].mean()
                    temp = temp[:,np.newaxis].T

                    # Shift template
                    temp = np.roll(temp, dt, 1)

                    # Same operation as to data
                    temp = temp - temp[0,-1]
                    
                    # HI templates
                    a = temp
                    s = 5
                    e = 120
                    if aa is None:
                        aa = a[:, s:e]
                        vv = v[s:e]
                        nsamp = len(vv)
                    else:
                        aa = np.hstack((aa, a[:, s:e]))
                        vv = np.hstack((vv, v[s:e]))


                # Add polynomial terms
                if mod is None:
                    mod = np.zeros((len(fwhmarr),len(darr),len(dtarr),len(dfarr),len(vv)))

                t = np.linspace(-1,1,nsamp)
                nterms = 2
                nb = len(barr)
                app = np.zeros((nterms*nb, len(vv)))
                for bb in range(nb):
                    s = bb*nsamp
                    e = (bb+1)*nsamp
                    ap = np.array([t**po for po in range(nterms)])
                    sb = bb*nterms
                    eb = (bb+1)*nterms
                    app[sb:eb,s:e] = ap
                aa = np.vstack((aa,app))

                # Linear least squares
                indd = np.isfinite(vv)
                x = np.linalg.lstsq(aa[:,indd].T, vv[indd])[0]

                if coef is None:
                    coef = np.zeros((len(fwhmarr),len(darr),len(dtarr),len(dfarr),len(x)))
                
                pred = aa.T.dot(x)

                mod[i,j,k,l,:] = pred
                chi2[i,j,k,l] = nansum((vv-mod[i,j,k,l,:])**2)
                coef[i,j,k,l,:] = x
                
ind = where(chi2 == np.min(chi2))
x = coef[ind]

fwhm = fwhmarr[ind[0][0]]
d = darr[ind[1][0]]
dt = dtarr[ind[2][0]]
df = dfarr[ind[3][0]]
print('fwhm={:f}, ddec={:f}, dt={:f}, df={:f}'.format(fwhm,d,dt,df))

clf()
plot(vv,label='data')
plot(mod[ind[0][0],ind[1][0],ind[2][0],ind[3][0],:],label='model')
xlabel('time index (1 idx = 1 deg RA)')
ylabel('T (K)')
ylim(-5,30)
legend(loc='lower right',prop={'size': 8})
for k,val in enumerate(barr):
    idx1 = (k)*nsamp
    idx2 = (k+1)*nsamp
    plot([idx1,idx1],[-5,30],':k')
    plot([idx2,idx2],[-5,30],':k')
    s = '{:0.2f}-{:0.2f}'.format(be_lo[val],be_hi[val])
    text((idx1+idx2)/2., 28, s, horizontalalignment='center',fontsize=8)


dra = (r.ra[2]-r.ra[1])*dt
title('fwhm={:0.1f}, ddec={:0.2f}, dra={:0.1f}, df={:0.2f}'.format(fwhm,d,dra,df))

#savefig('personal/hifitplots/hifit_{:s}.png'.format(r.tags[0][0:6]))

# Save best fit
#savez('personal/hifitplots/hifitcoefs_{:s}.npz'.format(r.tags[0][0:6]),
#      fwhm=fwhm, d=d, dt=dt, df=df, x=x)




