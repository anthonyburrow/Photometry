import os.path
import numpy as np


class Zeropoint:

    def __init__(self, date, aperture):
        self.date = date
        self.aperture = aperture

    def StandsList(self):
        root = "../../standards/"

        files = []
        if os.listdir(root) != []:
            if os.path.isdir(root + self.date):
                if os.listdir(root + self.date + "/") != []:
                    if os.path.isdir(root + self.date + "/ap" + str(self.aperture)):
                        if os.listdir(root + self.date + "/ap" + str(self.aperture) + "/") != []:
                            for file in sorted(os.listdir(root + self.date + "/ap" + str(self.aperture) + "/")):
                                if file[-6:] == ".mag.1":
                                    files.append(file)
        else:
            print("There is no standard star photometry available.")

        if files == []:
            print("No standard star photometry for " + self.date + " with aperture " + str(self.aperture))

        return files

    def magRead(self, filename):
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
                selected.extend((combined[7], combined[8], combined[33], combined[34]))
                data.append(selected)

        return data

    def ZeroPoint(self):
        filename = "../../standards/" + self.date + "/realMags.dat"
        # V, B-V, U-B, V-R, R-I
        realData = np.loadtxt(filename)

        files = self.StandsList()
        Bfiles = [file for file in files if file[0] == "B"]
        Vfiles = [file for file in files if file[0] == "V"]
        Rfiles = [file for file in files if file[0] == "R"]
        files = [Bfiles, Vfiles, Rfiles]

        # Make sure stars in .mag files are in the same order as realMag file here

        BSample = []
        VSample = []
        RSample = []

        for xFiles in files:
            # Set up filter-specific data set
            xData = []
            for file in xFiles:
                filename = "../../standards/" + self.date + "/ap" + self.aperture + "/" + file
                d = self.magRead(filename)
                if file == xFiles[0]:
                    numStars = len(d)
                xData.append(d)

            # For this filter, match all the same stars and add their zeropoint calculations to sample
            for i in range(0, numStars):
                star = []
                star.append(xData[i])
                for otherTarget in xData[numStars:]:
                    if abs(xData[i][0] - otherTarget[0]) <= 20 and abs(xData[i][1] - otherTarget[1]) <= 20:
                        star.append(otherTarget)

                for s in star:
                    if xFiles[0][0] == "B":
                        # zp_B = [ real(B-V) + real(V) ] - obs(B) + 25
                        zp = (realData[i][1] + realData[i][0]) - s[2] + 25
                        BSample.append(zp)
                    elif xFiles[0][0] == "V":
                        # zp_V = real(V) - obs(V) + 25
                        zp = realData[i][0] - s[2] + 25
                        VSample.append(zp)
                    elif xFiles[0][0] == "R":
                        # zp_R = [ real(V) - real(V-R) ] - obs(R) + 25
                        zp = (realData[i][0] - realData[i][3]) - s[2] + 25
                        RSample.append(zp)

        # Average sample to get mean zeropoint per filter
        BZeropoint = np.mean(BSample)
        VZeropoint = np.mean(VSample)
        RZeropoint = np.mean(RSample)

        # Write to file
        filename = "../../standards/" + self.date + "/ap" + self.aperture + "/zeropoints.dat"
        with open(filename, 'w') as F:
            data = [BZeropoint, VZeropoint, RZeropoint]
            np.savetxt(F, data, fmt='%.3f')
