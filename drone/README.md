Code for BMX beam mapping using drones
Note most code requires %pylab to load all from numpy and matplotlib.pyplot (laziness on my part)

droneDataPair.py	associates drone and BMX data, interpolates to BMX timestamps

fitOneDish.py		use scipy.optimize.curve_fit to fit 2D Gaussian to coordinate and signal data from droneDataPair()

Gaussfit_freqloop.py	example use of fitOneDish() in loop over dishes, frequecies

plottingCodeBMXmtg200504.py	example of how presentation slides were generated

visualizeBeam.py	example plot of interpolated beam overlaid with centroid from 2D Gaussian fit
