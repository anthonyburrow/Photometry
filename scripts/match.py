import numpy as np


class Match:
    """Creates a Match object, holding information on the star-matching process.

    Processes IRAF's output data into a more useable form and correlates (or matches)
    targets between different images for each filter and for each exposure time.  This
    outputs a data set describing the entire combined photometry for a single night
    of observation and star cluster.

    Attributes:
            output_directory: Directory desired for matched output.
            input_directory: Directory that holds IRAF's output photometry data.
            phot_type: Determines whether "psf" or "aperture" photometry is desired.
            coo_tolerance: Maximum pixel-distance acceptable for a match.
            mag_tolerance: Maximum magnitude difference acceptable for a match.
    """

    def __init__(self, output_directory, input_directory, phot_type="psf", coo_tolerance=5.0,
                 mag_tolerance=0.5):
        self.output_directory = output_directory
        self.input_directory = input_directory
        self.phot_type = phot_type
        self.coo_tol = coo_tolerance
        self.mag_tol = mag_tolerance

        self.short_psf_files = ["B1.als.1", "V1.als.1", "R1.als.1", "H1.als.1"]
        self.long_psf_files = ["B3.als.1", "V3.als.1", "R3.als.1", "H3.als.1"]
        self.short_aperture_files = ["B1.mag.1", "V1.mag.1", "R1.mag.1", "H1.mag.1"]
        self.long_aperture_files = ["B3.mag.1", "V3.mag.1", "R3.mag.1", "H3.mag.1"]

    def alsRead(self, filename):
        """Reads PSF photometry files.

        PSF photometry is in the standard output .als format provided by the 'allstar'
        task within IRAF's DAOPHOT package.  This file is located in
        root/photometry/*date*/*cluster*/

        Args:
                filename (string): The input .als file to read.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        F = open(self.input_directory + filename)
        file = F.readlines()[44:]
        F.close()

        data = []
        minimum = 0.
        maximum = 4096.
        for i in range(0, int(len(file) / 2.)):
            curr = file[2 * i].split()
            # Set image size limits (not needed if min/max equal image dimensions)
            if float(curr[1]) > minimum and float(curr[1]) < maximum and \
               float(curr[2]) > minimum and float(curr[2]) < maximum:
                # Concatenate lines
                combined = file[2 * i] + ' ' + file[2 * i + 1]
                # Select values needed in data set: X, Y, mag, mag error
                combined = combined.split()
                selected = []
                selected.extend((combined[1], combined[2], combined[3], combined[4]))
                data.append(selected)

        return data

    def magRead(self, filename):
        """Reads aperture photometry files.

        Aperture photometry is in the standard output .mag format provided by the 'phot'
        task within IRAF's DAOPHOT package.  This file is located in
        root/photometry/*date*/*cluster*/

        Args:
                filename (string): The input .mag file to read.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        F = open(self.input_directory + filename)
        file = F.readlines()[75:]
        F.close()

        data = []

        for i in range(0, int(len(file) / 5.)):
            # Concatenate lines
            combined = file[5 * i] + ' ' + file[5 * i + 1] + ' ' + file[5 * i + 2] + ' ' + file[5 * i + 3] + ' ' + file[5 * i + 4]
            # Select values needed in data set: X, Y, mag, mag error
            combined = combined.split()
            selected = []
            selected.extend((combined[7], combined[8], combined[33], combined[34]))
            data.append(selected)

        return data

    def ByFilter(self, exposure):
        """Matches targets between filters for each exposure time.

        The four filters used are B, V, R, and H-alpha.  Determines which targets have
        corresponding values within each filter data set.  Uses B filter data as a
        reference.

        Args:
                exposure (string): Determines whether "Short" or "Long" exposure times are used.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates and the
                magnitudes and magnitude errors for every target that has a corresponding
                value on each filter data set.

        """
        # Create data sets for each filter
        data = []

        if self.phot_type == "psf":
            if exposure == "Short":
                filenames = self.short_psf_files
            elif exposure == "Long":
                filenames = self.long_psf_files
            B_data = self.alsRead(filenames[0])
            V_data = self.alsRead(filenames[1])
            R_data = self.alsRead(filenames[2])
            H_data = self.alsRead(filenames[3])

        elif self.phot_type == "aperture":
            if exposure == "Short":
                filenames = self.short_aperture_files
            elif exposure == "Long":
                filenames = self.long_aperture_files
            B_data = self.magRead(filenames[0])
            V_data = self.magRead(filenames[1])
            R_data = self.magRead(filenames[2])
            H_data = self.magRead(filenames[3])

        # Specify any coordinate offsets left to be made (from the B image, which is the reference)
        V_coo_offset = [0, 0]
        R_coo_offset = [0, 0]
        H_coo_offset = [0, 0]

        # Match stars between filters
        for b in B_data:
            x_b = float(b[0])
            y_b = float(b[1])
            for v in V_data:
                x_v = float(v[0]) - V_coo_offset[0]
                y_v = float(v[1]) - V_coo_offset[1]
                if abs(x_b - x_v) < self.coo_tol and abs(y_b - y_v) < self.coo_tol:
                    for r in R_data:
                        x_r = float(r[0]) - R_coo_offset[0]
                        y_r = float(r[1]) - R_coo_offset[1]
                        if abs(x_b - x_r) < self.coo_tol and abs(y_b - y_r) < self.coo_tol:
                            for h in H_data:
                                x_h = float(h[0]) - H_coo_offset[0]
                                y_h = float(h[1]) - H_coo_offset[1]
                                if abs(x_b - x_h) < self.coo_tol and abs(y_b - y_h) < self.coo_tol:
                                    # Select values needed in data set: B_X, B_Y, B, Berr, V, Verr, R, Rerr, H, Herr
                                    selected = []
                                    selected.extend((b[0], b[1], b[2], b[3], v[2], v[3],
                                                     r[2], r[3], h[2], h[3]))
                                    data.append(selected)

        print("    " + exposure + " matched: " + str(len(data)))
        return data

    def ByExposure(self):
        """Matches targets between exposure times.

        Determines which targets have corresponding values between the short and
        long exposure data sets.  The magnitudes used for each filter are determined
        by those with the lesser respective error.  Targets with no matches are also
        included, as those typically refer to the brightest and darkest targets.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates and the
                magnitudes and magnitude errors of each filter for every target.

        """
        print("    Matching objects between filters...")

        # Create data sets for long and short exposures
        short_data = self.ByFilter("Short")
        long_data = self.ByFilter("Long")
        data = short_data + long_data

        print("    Matching objects between long and short exposures...")

        # Match between short and long exposures and use values from that with the lowest error
        count = 0
        for s in short_data:
            for l in long_data:
                if (abs(s[0] - l[0]) <= self.coo_tol) and (abs(s[1] - l[1]) <= self.coo_tol) and \
                   (abs(s[2] - l[2]) <= self.mag_tol) and (abs(s[4] - l[4]) <= self.mag_tol) and \
                   (abs(s[6] - l[6]) <= self.mag_tol) and (abs(s[8] - l[8]) <= self.mag_tol):

                    matched = [s[0], s[1]]
                    if (s[3] <= l[3]):
                        matched.extend((s[2], s[3]))
                    else:
                        matched.extend((l[2], l[3]))
                    if (s[5] <= l[5]):
                        matched.extend((s[4], s[5]))
                    else:
                        matched.extend((l[4], l[5]))
                    if (s[7] <= l[7]):
                        matched.extend((s[6], s[7]))
                    else:
                        matched.extend((l[6], l[7]))
                    if (s[9] <= l[9]):
                        matched.extend((s[8], s[9]))
                    else:
                        matched.extend((l[8], l[9]))

                    data.remove(s)
                    data.remove(l)
                    data.append(matched)

                    count += 1

        print
        print("    Matched between exposures: " + count)
        print("    Short only:                " + str(len(short_data) - count))
        print("    Long only:                 " + str(len(long_data) - count))
        print("    Total:                     " + str(len(data)))
        print

        # Output to file
        F = open(self.output_directory + "phot_" + self.phot_type + ".dat", 'w')

        for item in data:
            F.write(" ".join(item))
            F.write("\n")

        F.close()

        return data

    def LowError(self, max_error=0.04):
        """Determines which targets are within a margain of error.

        Determines the targets that exhibit a constrained error.  This outputs
        the full photometry for only targets within this error.

        Args:
                max_error: Maximum error that a R-H value is allowed.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates and the
                magnitudes and magnitude errors of each filter for every target that
                are within the given margain of error.

        """
        data = self.ByExposure()
        lowError_data = []

        Rerr = np.array(data[7])
        Herr = np.array(data[9])

        R_Herr = np.sqrt(Rerr**2 + Herr**2)

        for i in range(0, len(data)):
            if R_Herr[i] < max_error:
                lowError_data.append(data[i])

        # Output to file
        F = open(self.output_directory + "phot_" + self.phot_type + "_lowError.dat", 'w')

        for item in data:
            F.write(" ".join(item))
            F.write("\n")

        F.close()

        return lowError_data
