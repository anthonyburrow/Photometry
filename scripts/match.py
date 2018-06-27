import numpy as np


class Match:
    """Creates a Match object, holding information on the star-matching process.

    Processes IRAF's output data into a more useable form and correlates (or matches)
    targets between different images for each filter and for each exposure time.  This
    outputs a data set describing the entire combined photometry for a single night
    of observation and star cluster.

    Attributes:
            cluster:
            date:
            app:
    """

    def __init__(self, cluster, date, app):
        self.cluster = cluster
        self.date = date
        self.app = app

        self.short_psf_files = ["B1.als.1", "V1.als.1", "R1.als.1", "H1.als.1"]
        self.long_psf_files = ["B3.als.1", "V3.als.1", "R3.als.1", "H3.als.1"]
        self.short_aperture_files = ["B1.mag.1", "V1.mag.1", "R1.mag.1", "H1.mag.1"]
        self.long_aperture_files = ["B3.mag.1", "V3.mag.1", "R3.mag.1", "H3.mag.1"]

    def alsRead(self, file):
        """Reads PSF photometry files.

        PSF photometry is in the standard output .als format provided by the 'allstar'
        task within IRAF's DAOPHOT package.  This file is located in
        root/photometry/*date*/*cluster*/

        Args:
                file (string): The input .als file to read.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        filename = "../photometry/" + self.cluster + "/" + self.date + "/" + file
        with open(filename) as F:
            file = F.readlines()[44:]

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

    def magRead(self, file):
        """Reads aperture photometry files.

        Aperture photometry is in the standard output .mag format provided by the 'phot'
        task within IRAF's DAOPHOT package.  This file is located in
        root/photometry/*date*/*cluster*/

        Args:
                file (string): The input .mag file to read.

        Returns:
                2-dimensional array consisting of X- and Y- image coordinates, magnitudes,
                and magnitude errors.

        """
        filename = "../photometry/" + self.cluster + "/" + self.date + "/" + file
        with open(filename) as F:
            file = F.readlines()[75:]

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

        if self.app.phot_type == "psf":
            if exposure == "Short":
                filenames = self.short_psf_files
            elif exposure == "Long":
                filenames = self.long_psf_files
            B_data = self.alsRead(filenames[0])
            V_data = self.alsRead(filenames[1])
            R_data = self.alsRead(filenames[2])
            H_data = self.alsRead(filenames[3])

        elif self.app.phot_type == "aperture":
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
                if abs(x_b - x_v) < self.app.cooTol and abs(y_b - y_v) < self.app.cooTol:
                    for r in R_data:
                        x_r = float(r[0]) - R_coo_offset[0]
                        y_r = float(r[1]) - R_coo_offset[1]
                        if abs(x_b - x_r) < self.app.cooTol and abs(y_b - y_r) < self.app.cooTol:
                            for h in H_data:
                                x_h = float(h[0]) - H_coo_offset[0]
                                y_h = float(h[1]) - H_coo_offset[1]
                                if abs(x_b - x_h) < self.app.cooTol and abs(y_b - y_h) < self.app.cooTol:
                                    # Select values needed in data set: B_X, B_Y, B, Berr, V, Verr, R, Rerr, H, Herr
                                    selected = []
                                    selected.extend((float(b[0]), float(b[1]), float(b[2]), float(b[3]), float(v[2]), float(v[3]),
                                                     float(r[2]), float(r[3]), float(h[2]), float(h[3])))
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
        print("  Matching objects between filters...")

        # Create data sets for long and short exposures
        short_data = self.ByFilter("Short")
        long_data = self.ByFilter("Long")
        data = short_data + long_data

        print("  Matching objects between long and short exposures...")

        # Match between short and long exposures and use values from that with the lowest error
        count = 0
        for s in short_data:
            for l in long_data:
                if (abs(s[0] - l[0]) <= self.app.cooTol) and (abs(s[1] - l[1]) <= self.app.cooTol) and \
                   (abs(s[2] - l[2]) <= self.app.magTol) and (abs(s[4] - l[4]) <= self.app.magTol) and \
                   (abs(s[6] - l[6]) <= self.app.magTol) and (abs(s[8] - l[8]) <= self.app.magTol):

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

        print("    Matched between exposures: " + str(count))
        print("    Short only: " + str(len(short_data) - count))
        print("    Long only: " + str(len(long_data) - count))
        print("    Total: " + str(len(data)))

        # Output to file
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"
        with open(filename, 'w') as F:
            np.savetxt(F, data, fmt='%.3f')

        return data
