import numpy as np
from astrometry import Astrometry
from astropy import wcs
import csv


class Match:
    """Creates a Match object, holding information on the star-matching process.

    Processes IRAF's output data into a more useable form and correlates (or matches)
    targets between different images for each filter and for each exposure time.  This
    outputs a data set describing the entire combined photometry for a single night
    of observation and star cluster.

    Attributes:
            cluster: Cluster for which data to read from.
            date: Data for which data to read from.
            app: GUI application that specifies parameters.
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
        try:
            with open(filename) as F:
                file = F.readlines()[44:]
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

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
                selected.extend((float(combined[1]), float(combined[2]), float(combined[3]), float(combined[4])))
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
        try:
            with open(filename) as F:
                file = F.readlines()[75:]
        except IOError:
            print("\nFile does not exist:\n" + filename)
            return

        data = []

        for i in range(0, int(len(file) / 5.)):
            # Concatenate lines
            combined = file[5 * i] + ' ' + file[5 * i + 1] + ' ' + file[5 * i + 2] + ' ' + file[5 * i + 3] + ' ' + file[5 * i + 4]
            # Select values needed in data set: X, Y, mag, mag error
            combined = combined.split()
            selected = []
            if combined[33] != "INDEF" and combined[34] != "INDEF":
                selected.extend((float(combined[7]), float(combined[8]), float(combined[33]), float(combined[34])))
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
        if exposure == "Short":
            V_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="V1", baseImage="B1")
            R_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="R1", baseImage="B1")
            H_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="H1", baseImage="B1")
        if exposure == "Long":
            V_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="V3", baseImage="B3")
            R_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="R3", baseImage="B3")
            H_coo_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="H3", baseImage="B3")

        for target in V_data:
            target[0] += V_coo_offset[0]
            target[1] += V_coo_offset[1]
        for target in R_data:
            target[0] += R_coo_offset[0]
            target[1] += R_coo_offset[1]
        for target in H_data:
            target[0] += H_coo_offset[0]
            target[1] += H_coo_offset[1]

        # Match stars between filters
        data = self.FilterMatch(B_data, V_data)
        data = self.FilterMatch(data, R_data)
        data = self.FilterMatch(data, H_data)

        if exposure == "Short":
            for target in data:
                target.append("s")
        elif exposure == "Long":
            for target in data:
                target.append("l")

        print("\n    " + exposure + " matched: " + str(len(data)) + "\n")
        return data

    def FilterMatch(self, arr1, arr2):
        matches = []
        for i in arr1:
            potentialMatches = []
            for j in arr2:
                if abs(i[0] - j[0]) < self.app.cooTol and abs(i[1] - j[1]) < self.app.cooTol:
                    potentialMatches.append(j)

            # Match with the closest target
            if potentialMatches != []:
                d = [np.sqrt(i[0] - x[0])**2 + (i[1] - x[1])**2 for x in potentialMatches]
                d = np.array(d)
                match = potentialMatches[np.argmin(d)]
                matches.append(i + [match[2], match[3]])

        return matches

    def ExposureMatch(self, arr1, arr2):
        data = arr1 + arr2
        n_arr1 = len(arr1)
        n_arr2 = len(arr2)
        count = 0

        for i in arr1:
            potentialMatches = []
            for j in arr2:
                if abs(i[0] - j[0]) < self.app.cooTol and abs(i[1] - j[1]) < self.app.cooTol:
                    potentialMatches.append(j)

            # Match with the closest target
            if potentialMatches != []:
                d = [np.sqrt(i[0] - x[0])**2 + (i[1] - x[1])**2 for x in potentialMatches]
                d = np.array(d)
                match = potentialMatches[np.argmin(d)]

                exps = ""
                phot_used = [i[0], i[1]]

                for n in [0, 2, 4, 6]:
                    if (i[3 + n] >= match[3 + n]) and (abs(i[2 + n] - match[2 + n]) < self.app.magTol):
                        phot_used.extend((match[2 + n], match[3 + n]))
                        exps += "l"
                    else:
                        phot_used.extend((i[2 + n], i[3 + n]))
                        exps += "s"

                phot_used.append(exps)

                data.remove(i)
                data.remove(match)
                arr2.remove(match)  # Stop from matching again
                data.append(phot_used)

                count += 1

        print("\n    Matched between exposures: " + str(count))
        print("    Short only: " + str(n_arr1 - count))
        print("    Long only: " + str(n_arr2 - count))
        print("    Total: " + str(len(data)))

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
        print("  Matching objects between filters...\n")

        short_data = self.ByFilter("Short")
        long_data = self.ByFilter("Long")

        # Apply coordinate offset between B1 and B3
        coord_offset = Astrometry().GetOffset(self.cluster, self.date, baseDate=self.date, image="B3", baseImage="B1")
        for target in long_data:
            target[0] += coord_offset[0]
            target[1] += coord_offset[1]

        # Match between short and long exposures and use values from that with the lowest error
        print("\n  Matching objects between long and short exposures...")

        data = self.ExposureMatch(short_data, long_data)

        # Get RAs and DECs for distances
        w = wcs.WCS("../photometry/" + self.cluster + "/" + self.date + "/B1_wcs.fits")
        coords = []
        for target in data:
            radec = w.all_pix2world(target[0], target[1], 0)
            ra = float(radec[0])
            dec = float(radec[1])
            coords.append([ra, dec])

        filename = "../photometry/" + self.cluster + '/' + self.date + '/phot_radec.csv'
        with open(filename, 'w') as F:
            writer = csv.writer(F)
            writer.writerows(coords)

        # Correct for extinction:
        for target in data:
            target[2] -= self.app.A_b
            target[4] -= self.app.A_v
            target[6] -= self.app.A_r

        # Apply aperture correction
        try:
            filename = "../standards/" + self.date + "/" + self.cluster + "_aperture_corrections.dat"
            corrections = np.loadtxt(filename)
            for target in data:
                target[2] += corrections[0]
                target[4] += corrections[1]
                target[6] += corrections[2]
            print("  \nAperture corrections applied.")
        except IOError:
            print("  \nAperture corrections not applied.")

        # Output to file
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + ".dat"
        with open(filename, 'w') as F:
            np.savetxt(F, [x[:-1] for x in data], fmt='%.3f')

        filename = "../output/" + self.cluster + "/" + self.date + "/phot_specTypes_" + self.app.phot_type + ".dat"
        with open(filename, 'w') as F:
            for target in data:
                F.write('%8.3f' % target[0] + '    ')
                F.write('%8.3f' % target[1] + '    ')
                F.write('%6.3f' % target[2] + '    ')
                F.write('%6.3f' % target[4] + '    ')
                F.write('%6.3f' % target[6] + '    ')
                F.write('%6.3f' % target[8] + '    ')
                F.write(target[10])
                F.write("\n")

        return data
