# requires droneDataPair.py
# generates list of dicts, one per flight, with drone and BMX data
# compute-intensive; pickled results in 4.2GB file
# /gpfs01/astro/workarea/bmxdata/reduced/allFLYs.pkl

allFLYs = []
FLYS = ['340', '342', '343', '352', '353', '344', '347', '349', '351']
for i,f in enumerate(FLYS):
	allFLYs.append(droneDataPair(f, verbose=True, debug=True))
