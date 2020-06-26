# Plotting code for slide 4 in powerpoint
# requires %pylab, allFLYs data from 
# /gpfs01/astro/workarea/bmxdata/reduced/allFLYs.pkl
	    
#db43 = ['344','347','349','351']
db43 = [5,6,7,8]
#db30 = ['342','343','352','353']
db30 = [1,2,3,4]

# plotting code for slide 3 in powerpoint
clf()
spno=0
for i in range(4):		# FLY index
    for j in range(4):		# dish index
	spno += 1
	subplot(4,4,spno)
	print "i,j,subplot = ", i,j,spno
	fly=allFLYs[db43[i]]
	dx=degrees(fly['thetaX'][j])
	dy=degrees(fly['thetaY'][j])
	nshades=30
	z = log10(fly['autos'][j][:,50])
	tricontourf(dx,dy,z,nshades)
	tricontour(dx,dy,z,colors='k')
	grid(c='w', lw=1)
	axhline(0, c='k',lw=2)
	axvline(0, c='k',lw=2)
	axis('equal')
	xlim(-5,15)
	ylim(-10,10)
	title(fly['FLY']+" "+dishes[j]+"dish "+fly['pol']+" pol")
	suptitle("200312 flights with attn=30dB", fontsize='large')
	if (spno>12): xlabel('degrees')
	if ((spno % 4) == 1): ylabel('degrees')


# plotting code for slide 5
#######################
clf()
f0 = 1101.0742514133453    # cut0 starting frequency
df = 2.14843726            # cut0 frequency bin
for i,fr in enumerate([1110, 1200, 1250,1300, 1350,1400,1450,1500,1550]):
	subplot(3,3,i+1)
	findex = int((fr-f0)/df)
	fly=allFLYs[3]		# FLY353 30dB Epol Ntravel
	dx=degrees(fly['thetaX'][1])
	dy=degrees(fly['thetaY'][1])
	nshades=20
	#z = log10(fly['autos'][3][:,findex])
	z = fly['autos'][1][:,findex]
	#tricontourf(dx,dy,z,nshades, levels=[13,13.5,14,14.5,15,15.5,16,16.5])
	tricontourf(dx,dy,log10(z/max(z)),levels=arange(-3,0.1,.1))
	tricontour(dx,dy,z,colors='k')
	#colorbar()
	grid(c='w', lw=1)
	axhline(0, c='k',lw=2)
	axvline(0, c='k',lw=2)
	axis('equal')
	xlim(-5,15)
	ylim(-10,10)
	title(str(fr)+' MHz')
suptitle(fly['FLY']+' E dish E pol EW grid')

