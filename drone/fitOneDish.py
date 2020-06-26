from scipy.optimize import curve_fit

def twoD_Gaussian((x,y), amplitude, x0, y0, sigma_x, sigma_y, offset):
        """ 2D gaussian, 6 parameters parameters """
	x0=float(x0)
	y0=float(y0)
	g = offset + amplitude * np.exp(-( (x-x0)**2/2/sigma_x**2 + (y-y0)**2/2/sigma_y**2))
	return g.ravel()


def fitOneDish(dish, freq, FLYdata, verbose=False, debug=False):
        """ fit 2D Gaussian to DJI drone data
	    inputs::
	    dish	"Ndish", "Edish", "Sdish", "Wdish" or "Ndish2"
	    freq	frequency in MHz
	    FLYdata	dict with keys "x", "y", "z" (all in meters interpolated to mjdB
	    "thetaX", "thetaY" (radians), list, one row for each dish
	    "mjdD", "mjdB",
	    "autos" (list of 8 waterfalls)
	    "pol", "attn", "FLY", "BMXdir", run identifiers where "pol" = "N" or "E"
	 """

	pol = FLYdata['pol']
	mjdD = FLYdata['mjdD']
	mjdB = FLYdata['mjdB']
	thetaX = FLYdata['thetaX']
	thetaY = FLYdata['thetaY']
	#
	# ToDo: map pol to range of indices covering the portion of the flight in that polarization, only needed when polarization changed in flight. Maybe use 'yaw' from DJI todetermine polarization orientation

	# map dish to autospectrum channel:
	# NOTE: drondDataPair() selects only channels sensitive to transmitted pol
	if (dish == 'Ndish') or (dish == 'Ndish2'): chan=0
	if dish == 'Edish': chan=1
	if dish == 'Sdish': chan=2
	if dish == 'Wdish': chan=3
	    
        if (verbose):
	        print "Set chan = ", chan

	# get index to frequency array (for BMX "cut0" analysis)
	f0 = 1101.0742514133453    # cut0 starting frequency
	df = 2.14843726            # cut0 frequency bin
	findex = int((freq-f0)/df)
	if (verbose):
	        print "Set findex = ", findex
	    

	p0 = [3e16, 0, 0, 0.05, 0.05, 1e14]
	xdata = np.array([thetaX[chan], thetaY[chan]])
	ydata =   (FLYdata['autos'][chan][:,findex]) 

	bounds = ([1e12, -.15,-.15,0,0,1e10], [1e17,.25,.15,.1,.1,1e14])
	if (verbose):
	        print "starting fit for dish ",dish," pol ",pol, "freq ",freq
	po, pcov = curve_fit(twoD_Gaussian, xdata, ydata, p0=p0, bounds=bounds)

	if (verbose):
	        print "Done"
	if (debug):
	        return(po,pcov,xdata, ydata, twoD_Gaussian(xdata, *po))
        else:
		return ((po, pcov))

