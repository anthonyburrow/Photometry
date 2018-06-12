import numpy as np
from astropy.io import fits

class Scale:

	def __init__(self, cluster, phot_type = "psf", coo_tol = 10):
		self.cluster = cluster
		self.phot_type = phot_type
		self.coo_tol = coo_tol

		self.baseBin = 1
		self.baseData = BaseData()

	def Binning(date):
		with fits.open("../photometry/" + cluster + "/" + date + "/B1.fits") as file:
			bin = file[0].header["X_BINNING"] # Need to check

		return bin

	def Filter(data, filterTol = 20):

		filtered_data = []

	    for target in data:
	    	isAlone = true
	    	for otherTarget in [x for x in data if x != target]:
	    		r = np.sqrt((target[0] - otherTarget[0])**2 + (target[1] - otherTarget[1])**2)
	    		if r > filterTol:
	    			isAlone = false
	    			break
	    	if isAlone:
	    		filtered_data.append(target)

	    return filtered_data

	def BaseData():
		# Create base data and Filter (spacial filter)
		with open("../photometry/" + cluster + "/obs_dates.txt") as F:  # Establish first date as scaling base
			date = F.readline()
		data = np.loadtxt("../photometry/" + cluster + "/" + date + "/phot_" + phot_type + ".dat")
		data = Filter(data, 20)
		baseBin = Binning(date)

		# Remove outliers
		with open("../output/" + cluster + "/" + date + "/beList.dat") as filtered_data:
			for target in data:
				if target in filtered_data:
					data.remove(target)

		return data

	def Scale(date, xOffset, yOffset): # TODO: Call Scale and automate offsets
		# Create base data and Filter (spacial filter)
		orig_data = np.loadtxt("../photometry/" + cluster + "/" + date + "/phot_" + phot_type + ".dat")
		data = Filter(orig_data, 20)
		binning = Binning(date)

		# Remove outliers
		with open("../output/" + cluster + "/" + date + "/beList.dat") as filtered_data:
			for target in data:
				if target in filtered_data:
					data.remove(target)

		B_diff = []
		V_diff = []
		R_diff = []
		H_diff = []

		for target in data:
			for baseTarget in baseData:
				if abs(baseBin*baseTarget[0] - (binning*target[0] + xOffset)) <= coo_tol and
				   abs(baseBin*baseTarget[1] - (binning*target[1] + yOffset)) <= coo_tol:
					B_diff.append(baseTarget[2] - target[2])
					V_diff.append(baseTarget[4] - target[4])
					R_diff.append(baseTarget[6] - target[6])
					H_diff.append(baseTarget[8] - target[8])

		B_offset = np.mean(B_diff)
		V_offset = np.mean(V_diff)
		R_offset = np.mean(R_diff)
		H_offset = np.mean(H_diff)

		B_std = np.std(B_diff)
		V_std = np.std(V_diff)
		R_std = np.std(R_diff)
		H_std = np.std(H_diff)

		for target in orig_data:
			target[2] += B_offset
			target[4] += V_offset
			target[6] += R_offset
			target[8] += H_offset

			target[3] = np.sqrt(target[3]**2 + B_std**2)
			target[5] = np.sqrt(target[5]**2 + V_std**2)
			target[7] = np.sqrt(target[7]**2 + R_std**2)
			target[9] = np.sqrt(target[9]**2 + H_std**2)

		return orig_data