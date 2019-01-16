import numpy as np
from astropy.io import fits

import os.path


def GetAstrometryOffset(cluster, date, baseDate, image='B1', baseImage='B1'):
    """Calculates the coordinate offsets between two images.

    Images are specified by the dates and filter (file name) used for the image.

    Args:
        cluster (str): Cluster for the desired data.
        date (str): Variable date that is to be offset.
        baseDate (str): Reference date to be offset toward.
        image (str): Image name/filter of variable image
        baseImage (str): Image name/filter of reference image.

    Returns:
        list: List containing x- and y-offset to be applied.

    """
    # Get image values
    filename = 'photometry/' + cluster + '/' + baseDate + '/' + \
               baseImage + '.fits'
    F = fits.getheader(filename)
    baseMaxPixels = F['NAXIS1']
    baseBinning = F['XBINNING']

    filename = 'photometry/' + cluster + '/' + date + '/' + image + '.fits'
    F = fits.getheader(filename)
    binning = F['XBINNING']

    # Read plate scaled information
    baseFn = 'photometry/' + cluster + '/' + baseDate + '/' + \
             baseImage + '_corr.fits'
    fn = 'photometry/' + cluster + '/' + date + '/' + image + '_corr.fits'
    if not (os.path.isfile(baseFn) and os.path.isfile(fn)):
        # Ex.: The file
        # '/photometry/NGC663/20151102/cooOffsets_H3_to_B1_20151102.dat'
        # has the coord. offsets of H3 to B1 of the same data (only ones with
        # different dates are B1 to B1)
        filename = 'photometry/' + cluster + '/' + date + \
                   '/cooOffsets_' + image + '_to_' + baseImage + '_' + \
                   baseDate + '.dat'
        if not os.path.isfile(filename):
            print("\nCould not retrieve any coordinate offset information.  " +
                  "Using coordinate offset: [0, 0]")
            offsets = (0, 0)
            return offsets

        offsets = tuple(np.loadtxt(filename).tolist())
        print("\n  Using manual coordinate offset: \
              [%s, %s]" % offsets)
        return offsets

    with fits.open(baseFn) as baseFile, fits.open(fn) as file:
        baseCorr = baseFile[1].data
        corr = file[1].data

    # Calculate offsets
    count = 0
    sample = []
    for baseTarget in baseCorr:
        # Pick a star closer to the middle of the image to average out
        # exaggerated differences due to image rotation
        if baseMaxPixels * 0.1 <= baseTarget[0] <= baseMaxPixels * 0.9 and \
           baseMaxPixels * 0.1 <= baseTarget[1] <= baseMaxPixels * 0.9:
            for target in corr:
                # If they refer to the same index star:
                if baseTarget[6] == target[6] and baseTarget[7] == target[7]:
                    xOff = baseBinning * baseTarget[0] - binning * target[0]
                    yOff = baseBinning * baseTarget[1] - binning * target[1]
                    sample.append([xOff, yOff])
                    count += 1
        if count == 5:
            break

    offsets = tuple(np.mean(sample, axis=0).tolist())
    return offsets
