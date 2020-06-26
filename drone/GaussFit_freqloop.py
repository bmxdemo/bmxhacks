# requires %pylab, fitOneDish.py

def getTwoDfit(dish,pol,flower, fupper, npts, verbose=False):
	""" return 2D Gaussian fit with uncertainties over specified frequencies """

	fr=[]; poo=[]; pstdd=[]
        for freqq in linspace(flower, fupper, npts):
	        po,pcov = fitOneDish(dish,pol,freqq,FLY331,yd,verbose=verbose)
		fr.append(freqq)
		poo.append(po)
		pstdd.append(sqrt(diag(pcov)))
	return(array(fr), array(poo), array(pstdd))


# following code looks at all 4 x 2 dishes x polarizations across full frequency range and makes plots

"""
# run getTwoDfit()  for all dishes across 100 frequency points
f0,po0,pstd0 = getTwoDfit('Ndish','N',1110,1550,100, verbose=True)
f1,po1,pstd1 = getTwoDfit('Edish','N',1110,1550,100, verbose=True)
f2,po2,pstd2 = getTwoDfit('Sdish','N',1110,1550,100, verbose=True)
f3,po3,pstd3 = getTwoDfit('Wdish','N',1110,1550,100, verbose=True)
f4,po4,pstd4 = getTwoDfit('Ndish','E',1110,1550,100, verbose=True)
f5,po5,pstd5 = getTwoDfit('Edish','E',1110,1550,100, verbose=True)
f6,po6,pstd6 = getTwoDfit('Sdish','E',1110,1550,100, verbose=True)
f7,po7,pstd7 = getTwoDfit('Wdish','E',1110,1550,100, verbose=True)

# Get median pointing for each dish
mpox0=degrees(median(po0[5:75,1]))
mpoy0=degrees(median(po0[5:75,2]))
spox0=degrees(std(po0[5:75,3]))
spoy0=degrees(std(po0[5:75,4]))

mpox1=degrees(median(po1[5:75,1]))
mpoy1=degrees(median(po1[5:75,2]))
spox1=degrees(std(po1[5:75,3]))
spoy1=degrees(std(po1[5:75,4]))

mpox2=degrees(median(po2[5:75,1]))
mpoy2=degrees(median(po2[5:75,2]))
spox2=degrees(std(po2[5:75,3]))
spoy2=std(degrees(po2[5:75,4]))

mpox3=degrees(median(po3[5:75,1]))
mpoy3=degrees(median(po3[5:75,2]))
spox3=std(degrees(po3[5:75,3]))
spoy3=std(degrees(po3[5:75,4]))

mpox4=degrees(median(po4[5:75,1]))
mpoy4=degrees(median(po4[5:75,2]))
spox4=degrees(std(po4[5:75,3]))
spoy4=degrees(std(po4[5:75,4]))

mpox5=degrees(median(po5[5:75,1]))
mpoy5=degrees(median(po5[5:75,2]))
spox5=degrees(std(po5[5:75,3]))
spoy5=degrees(std(po5[5:75,4]))

mpox6=degrees(median(po6[5:75,1]))
mpoy6=degrees(median(po6[5:75,2]))
spox6=degrees(std(po6[5:75,3]))
spoy6=degrees(std(po6[5:75,4]))

mpox7=degrees(median(po7[5:75,1]))
mpoy7=degrees(median(po7[5:75,2]))
spox7=degrees(std(po7[5:75,3]))
spoy7=degrees(std(po7[5:75,4]))

print '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f' % (mpox0, mpox1, mpox2, mpox3, mpox4, mpox5, mpox6, mpox7); print '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f' % (mpoy0, mpoy1, mpoy2, mpoy3, mpoy4, mpoy5, mpoy6, mpoy7);

print '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f' % (spox0, spox1, spox2, spox3, spox4, spox5, spox6, spox7); print '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f' % (spoy0, spoy1, spoy2, spoy3, spoy4, spoy5, spoy6, spoy7);

# plotting code: beam widths vs frequency, all dishes 4x2 subplots
subplot(421)
plot(f0[:85], degrees(po0[:85,3]), '.')
plot(f0[:85], degrees(po0[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'NX', fontsize='small')

subplot(423)
plot(f1[:85], degrees(po1[:85,3]), '.')
plot(f1[:85], degrees(po1[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'EY', fontsize='small')

subplot(425)
plot(f2[:85], degrees(po2[:85,3]), '.')
plot(f2[:85], degrees(po2[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'SX', fontsize='small')

subplot(427)
plot(f3[:85], degrees(po3[:85,3]), '.')
plot(f3[:85], degrees(po3[:85,4]),'.'); grid(); ylim(1,3);text(1150,2.6,'WY', fontsize='small')

subplot(422)
plot(f4[:85], degrees(po4[:85,3]), '.')
plot(f4[:85], degrees(po4[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'NY', fontsize='small')

subplot(424)
plot(f5[:85], degrees(po5[:85,3]), '.')
plot(f5[:85], degrees(po5[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'EX', fontsize='small')

subplot(426)
plot(f6[:85], degrees(po6[:85,3]), '.')
plot(f6[:85], degrees(po6[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'SY', fontsize='small')

subplot(428)
plot(f7[:85], degrees(po7[:85,3]), '.')
plot(f7[:85], degrees(po7[:85,4]),'.'); grid(); ylim(1,3); text(1150,2.6,'WX', fontsize='small')

suptitle('sigmaX, sigmaY vs. freq; FLY331 200311_1841', fontsize='large')

* plotting code for pointing
>>> plot degrees(po0[:85,1]), degrees(po0[:85,2]), 'o', label='NX'
>>> plot degrees(po0[:75,1]), degrees(po0[:75,2]), 'o', label='NX'
>>> plot degrees(po1[:75,1]), degrees(po1[:75,2]), 'o', label='EY'
>>> plot degrees(po1[5:75,1]), degrees(po1[5:75,2]), 'o', label='EY'
>>> plot degrees(po0[5:75,1]), degrees(po0[5:75,2]), 'o', label='NX'
>>> plot degrees(po1[5:75,1]), degrees(po1[5:75,2]), 'o', label='EY'
>>> plot degrees(po2[5:75,1]), degrees(po2[5:75,2]), 'o', label='SX'
>>> plot degrees(po3[5:75,1]), degrees(po3[5:75,2]), 'o', label='WY'
>>> plot degrees(po4[5:75,1]), degrees(po4[5:75,2]), '^', label='NY'
>>> plot degrees(po5[5:75,1]), degrees(po5[5:75,2]), '^', label='EX'
>>> plot degrees(po6[5:75,1]), degrees(po6[5:75,2]), '^', label='SY'
>>> plot degrees(po7[5:75,1]), degrees(po7[5:75,2]), '^', label='WX
"""
