import numpy as np
from astropy.io import fits
from observations import Observations


class Astrometry:
    """Creates an AstrometryOffset object, holding information on calculating astrometry offsets.

    Uses astrometry.net to retrieve accurate spacial data for any given data.  This information
    is then used to calculate more precise offsets between observation nights.
    """

    def __init__(self):
        # For ARCSAT, pixel scale is 0.465"/pi
        self.pixel_scale = 0.465    # arcsec/pixel

    def GetOffset(self, cluster, date):
        """Calculates the coordinate offsets between dates.

        Uses the first night of observation for a cluster as the reference date.

        Args:
            cluster: Cluster for the desired data.
            date: Date of the observed data.

        """
        # Get reference date
        baseDate = Observations().ListDates(cluster)[0]

        with fits.open("../photometry/" + cluster + "/" + baseDate + "/B1.fits") as file:
            baseMaxPixels = file[0].header['NAXIS1']
            baseBinning = file[0].header['XBINNING']
        with fits.open("../photometry/" + cluster + "/" + date + "/B1.fits") as file:
            binning = file[0].header['XBINNING']

        # Read new information
        try:
            with fits.open("../photometry/" + cluster + "/" + baseDate + "/B1_corr.fits") as file:
                baseCorr = file[1].data
        except IOError:
            print("\nRetrieve 'B1_corr.fits' file for " + cluster + " on " + baseDate + " before calculating astrometry offsets.")
            offsets = [0, 0]
            return offsets

        try:
            with fits.open("../photometry/" + cluster + "/" + date + "/B1_corr.fits") as file:
                corr = file[1].data
        except IOError:
            print("\nRetrieve 'B1_corr.fits' file for " + cluster + " on " + date + " before calculating astrometry offsets.")
            offsets = [0, 0]
            return offsets

        # Get accurate RA and DEC values
        # with fits.open("../photometry/" + cluster + "/" + baseDate + "/") as file:
        #     baseRA = file[0].header["CRVAL1"]
        #     baseDec = file[0].header["CRVAL2"]
        #
        # with fits.open("../photometry/" + cluster + "/" + date + "/") as file:
        #     RA = file[0].header["CRVAL1"]
        #     Dec = file[0].header["CRVAL2"]

        # Calculate offsets
        count = 0
        sample = []
        for baseTarget in baseCorr:
            if baseTarget[0] >= baseMaxPixels * 0.25 and baseTarget[0] <= baseMaxPixels * 0.75 and \
               baseTarget[1] >= baseMaxPixels * 0.25 and baseTarget[1] <= baseMaxPixels * 0.75:
                for target in corr:
                    if baseTarget[6] == target[6] and baseTarget[7] == target[7]:
                        xOff = baseBinning * baseTarget[0] - binning * target[0]
                        yOff = baseBinning * baseTarget[1] - binning * target[1]
                        sample.append([xOff, yOff])
                        count += 1
            if count == 5:
                break

        sample = np.array(sample)
        offsets = [np.mean(sample[:, 0]), np.mean(sample[:, 1])]
        offsets = np.array(offsets)
        return offsets
