import numpy as np
from astropy.io import fits


class Astrometry:
    """Creates an AstrometryOffset object, holding information on calculating astrometry offsets.

    Uses astrometry.net to retrieve accurate spacial data for any given data.  This information
    is then used to calculate more precise offsets between observation nights.
    """

    def __init__(self):
        pass

    def GetOffset(self, cluster, date, baseDate, image="B1", baseImage="B1"):
        """Calculates the coordinate offsets between dates.

        Uses the first night of observation for a cluster as the reference date.

        Args:
            cluster: Cluster for the desired data.
            date: Date of the observed data.
            baseDate:
            image:
            baseImage:

        """
        print("  Calculating coordinate offsets of " + date + ", " + image + ".fits from " + baseDate + ", " + baseImage + ".fits...")

        # Get image values
        F = fits.getheader("../photometry/" + cluster + "/" + baseDate + "/" + image + ".fits")
        baseMaxPixels = F['NAXIS1']
        baseBinning = F['XBINNING']

        F = fits.getheader("../photometry/" + cluster + "/" + date + "/" + image + ".fits")
        binning = F['XBINNING']

        # Read plate scaled information
        try:
            with fits.open("../photometry/" + cluster + "/" + baseDate + "/" + baseImage + "_corr.fits") as file:
                baseCorr = file[1].data
        except IOError:
            print("\nError: Retrieve '" + baseImage + "_corr.fits' file for " + cluster + " on " + baseDate + " before calculating astrometry offsets.")
            offsets = [0, 0]
            return offsets

        try:
            with fits.open("../photometry/" + cluster + "/" + date + "/" + image + "_corr.fits") as file:
                corr = file[1].data
        except IOError:
            print("\nError: Retrieve '" + image + "_corr.fits' file for " + cluster + " on " + date + " before calculating astrometry offsets.")
            offsets = [0, 0]
            return offsets

        # Calculate offsets
        count = 0
        sample = []
        for baseTarget in baseCorr:
            # Pick a star closer to the middle of the image to average out exaggerated differences due to image rotation
            if baseTarget[0] >= baseMaxPixels * 0.25 and baseTarget[0] <= baseMaxPixels * 0.75 and \
               baseTarget[1] >= baseMaxPixels * 0.25 and baseTarget[1] <= baseMaxPixels * 0.75:
                for target in corr:
                    if baseTarget[6] == target[6] and baseTarget[7] == target[7]:   # If they refer to the same index star
                        xOff = baseBinning * baseTarget[0] - binning * target[0]
                        yOff = baseBinning * baseTarget[1] - binning * target[1]
                        sample.append([xOff, yOff])
                        count += 1
            if count == 5:
                break

        sample = np.array(sample)
        offsets = [np.mean(sample[:, 0]), np.mean(sample[:, 1])]
        print("    Offset found to be ", offsets)
        offsets = np.array(offsets)
        return offsets
