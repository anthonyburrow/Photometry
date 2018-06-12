import numpy as np

class BeFilter:
    """Creates a BeFilter object, holding information on determining Be candidates.

    Be candidates are extracted and output into a single respective file.  The threshold
    that determines the candidates can be calculated automatically or specified manually.

    Attributes:
        output_directory: Directory desired for matched output.
        phot_type: Determines whether "psf" or "aperture" photometry is desired.
        auto_limit: When true, the R-H limit will be calculated automatically.
        R_H_threshold: Manually issue a specific threshold for determining candidates.
        B_V_limits: Specifies the color range desired for the filter process.
    """

	def __init__(self, output_directory, phot_type = "psf", auto_limit = true, 
				 R_H_threshold = -3.75, B_V_limits = [-0.5, 1]):
		self.phot_type = phot_type
		self.R_H_threshold = R_H_threshold
		self.B_V_limits = B_V_limits
		self.output_directory = output_directory

	def Full():
    """Extracts Be candidate data for every target.

    Issues the command to filter the entirety of the finalized photometry data.
    """
		data = np.loadtxt(output_directory + "phot_" + phot_type + ".dat")

		Filter(data, "beList.dat")

	def LowError():
    """Extracts Be candidate data for targets with lower error.

    Issues the command to filter the finalized photometry data which exhibits
    constrained error.
    """
		data = np.loadtxt(output_directory + "phot_" + phot_type + "_lowError.dat")

		Filter(data, "beList_lowError.dat")

	def Filter(data, output):
        """Determines which targets lie outside the threshold.

        Determines which targets lie outside the threshold and writes to a corresponding
        output.

        Args:
            data (array): Data set that is to be filtered.
            output (string): The filename for the output data.

        Returns:
            2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
            and magnitude errors for each target that is filtered.

        """
		B = data[:,2]
		V = data[:,4]
		R = data[:,6]
		H = data[:,8]

		B_V = B - V			# may not work, may need to make numpy arrays explicitly
		R_H = R - H

		if auto_limit:
			R_H_threshold = AutoThreshold(R_H)

		filtered_data = []
		for i in range(0, len(data)):
			if R_H[i] > R_H_threshold and B_V[i] > B_V_limits[0] and B_V[i] < B_V_limits[1]:
				filtered_data.append(data[i])
	    
	    # Output to file
	    F = open(output_directory + output,'w')

	    for item in filtered_data:
	    	F.write(" ".join(item))
	    	F.write("\n")

	    F.close()

	    return filtered_data

	def AutoThreshold(r_h, iterate_limit = 10):
		"""Automatically determines the R-H threshold.

        Statisically calculates the R-H threshold by iteratively deciding which targets
        lie outside three-sigma (3 times the standard deviation) of the mean R-H value.

        Args:
            r_h (array): R-H data for which the threshold is desired.
            iterate_limit: Maximum number of iterations that is allowed to calculate the limit.

        Returns:
            The newly calculated threshold value.

        """
        print ("	Calculating threshold...")

		mean = np.mean(r_h)
		std = np.std(r_h)
		threshold = mean + 3*std

		count = 0
		while count < iterate_limit:
			new_r_h = []
			for value in r_h:
				if value <= (mean + 3*std) and value >= (mean - 3*std):
					new_r_h.append(value)

			new_mean = np.mean(new_r_h)
			new_std = np.std(new_r_h)
			threshold = new_mean + 3*new_std

			if mean != new_mean and std != new_std: 
				r_h = new_r_h
				mean = new_mean
				std = new_std
			else:
				break

			count += 1

		return threshold		