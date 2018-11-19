from astropy.io import fits
from astropy import wcs

import warnings


def Binning(cluster, date):
    """Determines the bin size of an image.

    Reads the B1.fits image for the specified date of observation and extracts
    the binning data.  This uses x-axis binning, however x- and y- axis binning
    are typically equal.

    Args:
            cluster (str): Cluster from which the image is retrieved.
            date (str): Date from which the image is retrieved.

    Returns:
            int: The bin size of the image.

    """
    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    try:
        F = fits.getheader(filename)
        binning = F['XBINNING']
    except IOError:
        print("\nFile does not exist:\n" + filename)
        binning = 1

    return binning


def GetFITSValues(cluster, date, keywords):
    result = {}

    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    try:
        for k in keywords:
            F = fits.getheader(filename)
            result[k] = F[k]
    except IOError:
        print("\nFile does not exist:\n" + filename)

    return result


def GetWCS(cluster, date):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        filename = 'photometry/' + cluster + '/' + date + \
                   '/B1_wcs.fits'
        try:
            w = wcs.WCS(filename)
        except IOError:
            print("\nError: Retrieve 'B1_wcs.fits' file for " +
                  cluster + " on " + date +
                  " before calculating exact RA and Dec values. \
                  Image center values added as placeholder.")

    return w


def alsRead(filename):
    """Reads PSF photometry files.

    PSF photometry is in the standard output .als format provided by the
    'allstar' task within IRAF's DAOPHOT package.  This file is located in
    `root/photometry/*date*/*cluster*/`.

    Args:
        filename (str): The input .als file (full path) to read.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates,
              magnitudes, and magnitude errors from the given file.

    """
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
        if minimum < float(curr[1]) < maximum and \
           minimum < float(curr[2]) < maximum:
            combined = ' '.join([file[2 * i], file[2 * i + 1]])
            # Select values needed in data set: X, Y, mag, mag error
            combined = combined.split()
            selected = []
            selected.extend([float(combined[1]), float(combined[2]),
                             float(combined[3]), float(combined[4])])
            data.append(selected)

    return data


def magRead(filename):
    """Reads aperture photometry files.

    Aperture photometry is in the standard output .mag format provided by the
    phot' task within IRAF's DAOPHOT package.  This file is located in
    root/photometry/*date*/*cluster*/

    Args:
        filename (str): The input .mag file (full path) to read.

    Returns:
        list: 2-dimensional array consisting of x- and y-image coordinates,
              magnitudes, and magnitude errors.

    """
    try:
        with open(filename) as F:
            file = F.readlines()[75:]
    except IOError:
        print("\nFile does not exist:\n" + filename)
        return

    data = []

    for i in range(0, int(len(file) / 5.)):
        # Concatenate lines
        combined = ' '.join(file[5 * i: 5 * i + 5])
        # Select values needed in data set: X, Y, mag, mag error
        combined = combined.split()
        selected = []
        if combined[33] != 'INDEF' and combined[34] != 'INDEF':
            selected.extend((float(combined[7]), float(combined[8]),
                             float(combined[33]), float(combined[34])))
            data.append(selected)

    return data
