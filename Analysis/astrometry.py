import numpy as np
from astropy.io import fits

import os.path


def GetAstrometryOffset(cluster, date, baseDate, image='B1', baseImage='B1'):
    """Calculates the coordinate offsets between dates.

    Uses the first night of observation for a cluster as the reference date.

    Args:
        cluster: Cluster for the desired data.
        date: Date of the observed data.
        baseDate:
        image:
        baseImage:

    """
    # Get image values
    filename = 'photometry/' + cluster + '/' + baseDate + '/' + baseImage + '.fits'
    F = fits.getheader(filename)
    baseMaxPixels = F['NAXIS1']
    baseBinning = F['XBINNING']

    filename = 'photometry/' + cluster + '/' + date + '/' + image + '.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

    # Read plate scaled information
    baseFn = 'photometry/' + cluster + '/' + baseDate + '/' + baseImage + '_corr.fits'
    fn = 'photometry/' + cluster + '/' + date + '/' + image + '_corr.fits'
    if not (os.path.isfile(baseFn) and os.path.isfile(fn)):
        # Ex.: The file '/photometry/NGC663/20151102/cooOffsets_H3_to_B1_20151102.dat'
        # has the coord. offsets of H3 to B1 of the same data (only ones with
        # different dates are B1 to B1)
        filename = 'photometry/' + cluster + '/' + date + '/cooOffsets_' + image + '_to_' + baseImage + '_' + baseDate + '.dat'
        if not os.path.isfile(filename):
            print("\nCould not retrieve any coordinate offset information.  " +
                  "Using coordinate offset: [0, 0]")
            offsets = [0, 0]
            return offsets

        offsets = np.loadtxt(filename).tolist()
        print("\n  Using manual coordinate offset: [%s, %s]" % (offsets[0], offsets[1]))
        return offsets

    with fits.open(baseFn) as baseFile, fits.open(fn) as file:
        baseCorr = baseFile[1].data
        corr = file[1].data

    # Calculate offsets
    count = 0
    sample = []
    for baseTarget in baseCorr:
        # Pick a star closer to the middle of the image to average out exaggerated differences due to image rotation
        if baseMaxPixels * 0.25 <= baseTarget[0] <= baseMaxPixels * 0.75 and \
           baseMaxPixels * 0.25 <= baseTarget[1] <= baseMaxPixels * 0.75:
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
    # print("    Offset found to be ", offsets)
    offsets = np.array(offsets)
    return offsets
