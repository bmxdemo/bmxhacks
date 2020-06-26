# interpolate mjd-Drone to mjd-BMX and calculate angles (W dish)
droneData = FLY349dd
BMXData = BMX2011
dish = "Ndish"
mjdD = droneData['mjd']
mjdB = BMXData['mjd']
yp = droneData['y']             # WAS yI
xp = droneData['x']
hmslp = droneData['hmsl']
y_interp = np.interp(mjdB, mjdD, yp)
x_interp = np.interp(mjdB, mjdD, xp)
hmsl_interp = np.interp(mjdB, mjdD, hmslp)
thetaX = np.arctan((x_interp - dishLoc[dish][1])/(hmsl_interp - dishLoc[dish][2]))
thetaY = np.arctan((y_interp - dishLoc[dish][0])/(hmsl_interp - dishLoc[dish][2]))
#
########
# or 3D surface plot
from mpl_toolkits import mplot3d
norm = matplotlib.colors.Normalize(vmin = np.min(bigz), vmax = np.max(bigz), clip = False)
ax = plt.axes(projection='3d')
ax.plot_surface(xx,yy,(tryOut), cmap='rainbow',norm=norm, rstride=4, cstride=4, edgecolor='none')


#######################
for d in ["Ndish", "Edish", "Sdish", "Wdish"]:
	for f in allFLYs:
		po, pcov = fitOneDish(d, 1350, f, dishLoc)
		po3 = array_str(degrees(po[1:5]),precision=3)
		pc3 = array_str(degrees(sqrt(diag(pcov))[1:5]),precision=3)
		print f["FLY"], f["pol"], d, po3,pc3

#CREATE FIGURE slide 2 BMXmgt_200504.pptx
# 2D beam linear scale using griddata
# generate grid of points covering theta 
clf()
mgslice = slice(15000,40000)   # abbreviated set of xy points
npts = 200
xx,yy = meshgrid(linspace(-.20,.20,npts), linspace(-.20,.20,npts))
points = (xdata[0][mgslice], xdata[1][mgslice])
values = ydata[mgslice]
tryOut = griddata(points, values, (xx,yy), rescale=True, method='nearest')
imshow(tryOut, origin='lower', interpolation='nearest')

####################

dishLoc = np.array([[  4.53866934e+06,  -6.12354322e+06,   1.70000000e+01],
[  4.53866568e+06,  -6.12353927e+06,   1.70000000e+01],
[  4.53866241e+06,  -6.12354284e+06,   1.70000000e+01],
[  4.53866555e+06,  -6.12354663e+06,   1.70000000e+01]])
######################
r351 = r'https://docs.google.com/spreadsheets/d/e/2PACX-1vTplEvCskI_GO3yb7YDNFm8thTfK4mmoYpSh-s8anNRqQ5Yjg2ZLCYvEKl22uwz1FHTd60M2Lx9vDrR/pub?output=csv'
FLY351=pd.read_csv(r3351)
FLY351d = FLY351.to_dict(orient = 'list')

##########################
t = Time(inDict['datetimestamp'])
    MJD0 = t.mjd[0]
    se = array(inDict['sec_elapsed']) - inDict['sec_elapsed'][0]
    mjd = MJD0 + se/86400.
##########################
clf()
for i in range(4):
    subplot(2,2,i+1)
    dx=degrees(jnk['thetaX'][i])
    dy=degrees(jnk['thetaY'][i])
    nshades=20
    z = log10(jnk['autos'][i][:,110])
    tricontourf(dx,dy,z,nshades)
    tricontour(dx,dy,z,colors='k')
    grid(c='w')
    axhline(0, c='k')
    axvline(0, c='k')
    xlim(-5,15)
    ylim(-10,10)
    axis('equal')

###########################
# generate list of dicts for all FLYs
# pickled in ...workarea/bmxdata/reduced/fra

DJI200312runs = ['340','342','343','352','353','344','347','349','351']
allFLYs=[]
for FLY in DJI200312runs:
    allFLYs.append(droneDataPair(FLY, verbose=True))


###########################
# CREATE FIGURE slide 3 + 4 BMXmgt_200504.pptx
# 2D interp log beam 4 FLYs, 4 dishes
dishes = ['N','E','S','W']
#db43 = ['344','347','349','351']
db43 = [5,6,7,8]
#db30 = ['342','343','352','353']
db30 = [1,2,3,4]

clf()
spno=0
for i in range(4):		# FLY index
    for j in range(4):		# dish index
	spno += 1
	subplot(4,4,spno)
	print "i,j,subplot = ", i,j,spno
	#fly=allFLYs[db30[i]]
	fly=allFLYs[db43[i]]
	dx=degrees(fly['thetaX'][j])
	dy=degrees(fly['thetaY'][j])
	nshades=64
	z = log10(fly['autos'][j][:,110])
	tricontourf(dx,dy,z,nshades)
	tricontour(dx,dy,z,colors='k')
	grid(c='w', lw=1)
	axhline(0, c='k',lw=2)
	axvline(0, c='k',lw=2)
	axis('equal')
	xlim(-5,15)
	ylim(-10,10)
	title(fly['FLY']+" "+dishes[j]+"dish "+fly['pol']+" pol")
	suptitle("200312 flights with attn= "+str(fly['attn'])+"dB", fontsize='large')
	if (spno>12): xlabel('degrees')
	if ((spno % 4) == 1): ylabel('degrees')

###########################

#CREATE FIGURE slide 5 BMXmgt_200504.pptx
# 2D interp log beam vs. freq, one dish, onr FLY
clf()
f0 = 1101.0742514133453    # cut0 starting frequency
df = 2.14843726            # cut0 frequency bin
for i,fr in enumerate([1110, 1200, 1250,1300, 1350,1400,1450,1500,1550]):
	subplot(3,3,i+1)
	findex = int((fr-f0)/df)
	# but looks like figure was produced for FLY352, allFLYs[8]
	fly=allFLYs[5]		# FLY344 43dB Npol Etravel
	dx=degrees(fly['thetaX'][3])
	dy=degrees(fly['thetaY'][3])
	#nshades=20
	#z = log10(fly['autos'][3][:,findex])
	z = fly['autos'][3][:,findex]
	#tricontourf(dx,dy,z,nshades, levels=[13,13.5,14,14.5,15,15.5,16,16.5])
	tricontourf(dx,dy,log10(z/max(z)),levels=arange(-4,0.1,.1))
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

#########################

# try concatenating EW- and NX-travel directions same pol, dish, freq
# 43dB files 344/347 (N pol, EW/NS <->) allFLYs 5,6
# 349/351 (E pol, EW/NS <->) allFLYs 7,8
# 30dB files 342/343 (N pol EW/NS <->), allFLYs 1,2

fly1 = allFLYs[7]; fly2=allFLYs[8] #349/351
#fly1 = allFLYs[1]; fly2=allFLYs[8] #349/351
dishindex=2
dx1 = fly1['thetaX'][dishindex]		# W dish
dx2 = fly2['thetaX'][dishindex]		# W dish
dy1 = fly1['thetaY'][dishindex]		# W dish
dy2 = fly2['thetaY'][dishindex]		# W dish
z1 = fly1['autos'][dishindex][:,50]	# W dish, 1208MHz
z2 = fly2['autos'][dishindex][:,50]	# W dish, 1208MHz
bigx=append(dx1,dx2)
bigy=append(dy1,dy2)
bigz=append(z1,z2)
#tricontour(bigx,bigy,log10(bigz), colors='k')
tricontourf(bigx,bigy,log10(bigz))

#########################
#CREATE FIGURE BMXmtg200511.pptx compare interp and 2DGauss fit

clf()

# index FLY POL DIR ATTN
#fly = allFLYs[1] #342  N E 30
#fly = allFLYs[2] #343  N N 30
fly = allFLYs[3] #352  E E 30
#fly = allFLYs[4] #353  E N 30
#fly = allFLYs[5] #344  N E 43
#fly = allFLYs[6] #347  N N 43
#fly = allFLYs[7] #349  E E 43
#fly = allFLYs[8] #351  E N 43

dishindex=0 # 0 N 1 E 2 S 3 W
hifly = fly['z'] > 175
dx = fly['thetaX'][dishindex][hifly]
dy = fly['thetaY'][dishindex][hifly]
z = fly['autos'][dishindex][:,50][hifly]
tricontourf(dx,dy,log10(z),10)
axis('equal')
xlim(-.20,.20); ylim(-.20,.20)
axhline(0, c='k')
axvline(0, c='k')
grid()

# To compare interpolated and 2D Gaussian fits
 po, pcov, xdata, ydata, fittedData  = fitOneDish(dishes[dishindex]+"dish", 1208, fly, dishLoc, debug=True)
plot ([po[1], po[1]], [po[2], po[2]], '^', markersize=12, c='k', label='Gaussfit center '+array_str(po[1:3], precision=3))
legend(loc='upper left', fontsize='small')
xlabel('radians')
ylabel('radians')
title(fly['FLY'] + " " + dishes[dishindex] + " dish " + fly['pol'] + "pol " + str(fly['attn']) + "dB attn 1208MHz")
# CREATE FIGURE
# beam for 4 dishes for one run at frequency index 50 (1208MHz)
# interpolated beam (contours) and center from 2D Gauss fit (triangles)
clf()
for dishindex in range(4):
    subplot(1,4,dishindex+1)
    hifly = fly['z'] > 175	# precaution to avoid start and end of flights
    dx = fly['thetaX'][dishindex][hifly]
    dy = fly['thetaY'][dishindex][hifly]
    z = fly['autos'][dishindex][:,50][hifly]
    print "Doing tricontourf for ", dishindex
    tricontourf(dx,dy,log10(z),10)
    axis('equal')
    xlim(-.20,.20); ylim(-.20,.20)
    axhline(0, c='k')
    axvline(0, c='k')
    grid()
    print "Doing fitOneDish for ", dishindex
    po, pcov, xdata, ydata, fittedData  = fitOneDish(dishes[dishindex]+"dish", 1208, fly, dishLoc, debug=True)
    plot ([po[1], po[1]], [po[2], po[2]], '^', markersize=12, c='k', label='Gaussfit center '+array_str(po[1:3], precision=3))
    legend(loc='upper left', fontsize='small')
    xlabel('radians')
    ylabel('radians')
    title(fly['FLY'] + " " + dishes[dishindex] + " dish " + fly['pol'] + "pol " + str(fly['attn']) + "dB attn 1208MHz")
    




#############################
#############################
# CREATE FIGURE
# get radial profile from griddata interpolate
clf()
hifly = fly['z'] > 175
dx = fly['thetaX'][dishindex][hifly]
dy = fly['thetaY'][dishindex][hifly]
z = fly['autos'][dishindex][:,50][hifly]
# with scipy.interpolate griddata (note matplotlib.mlab fcn same name)
clf()
npts=512
xx,yy = meshgrid(linspace(-.2, .2, npts), linspace(-.2,.2,npts))
points = (dx, dy)
tryOut = griddata(points, z, (xx,yy))
subplot(121)
#imshow(tryOut, interpolation='nearest', origin='lower')
imshow(log10(tryOut), interpolation='nearest', extent=[-.2,.2,-.2,.2], origin='lower')
grid()
xlabel('radians')
ylabel('radians')
title(fly['FLY']+" " + dishes[dishindex]+' dish ' + fly['pol'] + ' pol ' + str(fly['attn']) + ' dB')
subplot(122)
xprof = nanmean(tryOut,0)
xmax = argmax(yprof)
yprof = nanmean(tryOut,1)
ymax = argmax(yprof)
plot(xx[0], 10*log10(xprof), label='E-W')
plot(xx[0], 10*log10(yprof), label='N-S')
#plot(xx[0], xprof, label='E-W')
#plot(xx[0], yprof, label='N-S')
xlabel('radians')
ylabel('dB')
legend(loc='upper left')
grid()



