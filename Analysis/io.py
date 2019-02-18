from astropy.io import fits
from astropy import wcs

import warnings

import csv
import os.path


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
        print("\nFile does not exist:\n%s" % filename)
        binning = 1

    return binning


def GetFITSValues(cluster, date, keywords):
    """Gets specific values from a FITS header.

    Args:
        cluster (str): Cluster corresponding to the FITS image.
        date (str): Date corresponding to the FITS image.
        keywords (list): List of strings pertaining to header keywords desired.

    Returns:
        dictionary: Contains keys and values for the FITS information.

    """
    result = {}

    filename = 'photometry/' + cluster + '/' + date + '/B1.fits'
    try:
        for k in keywords:
            F = fits.getheader(filename)
            result[k] = F[k]
    except IOError:
        print("\nFile does not exist:\n%s" % filename)

    return result


def GetWCS(cluster, date):
    """Gets WCS image for a given observation.

    Args:
        cluster (str): Cluster corresponding to the WCS image desired.
        date (str): Date corresponding to the WCS image desired.

    Returns:
        wcs: The 'B1_wcs.fits' WCS image from a specific observation date.

    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        filename = 'photometry/' + cluster + '/' + date + \
                   '/B1_wcs.fits'
        try:
            w = wcs.WCS(filename)
        except IOError:
            print("\nError: Retrieve 'B1_wcs.fits' file for %s on %s before \
                   calculating exact RA and Dec values. Image center values \
                   added as placeholder." % (cluster, date))

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
        print("\nFile does not exist:\n%s" % filename)
        return

    n = 2   # Lines per target
    data = []
    for i in range(int(len(file) / n)):
        combined = ' '.join([file[n * i: n * (i + 1)]])
        combined = combined.split()

        # Select values needed in data set: X, Y, mag, mag error
        selected = []
        indices = (1, 2, 3, 4)
        for i in indices:
            selected.append(float(combined[i]))
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
        print("\nFile does not exist:\n%s" % filename)
        return

    n = 5   # Lines per target
    data = []
    for i in range(int(len(file) / n)):
        combined = ' '.join(file[n * i: n * (i + 1)])
        combined = combined.split()

        # Select values needed in data set: X, Y, mag, mag error
        selected = []
        indices = (7, 8, 33, 34)
        if combined[33] != 'INDEF' and combined[34] != 'INDEF':
            for i in indices:
                selected.append(float(combined[i]))
            data.append(selected)

    return data


def ReadClusterProperty(cluster, keys):
    filename = 'output/%s/cluster_properties.csv' % cluster
    reader = {}
    if os.path.isfile(filename):
        with open(filename, 'r') as F:
            reader = csv.DictReader(F)
            reader = [dict(row) for row in reader][0]

    values = []
    for key in keys:
        values.append(reader[key])

    return values


def WriteClusterProperty(cluster, key, value):
    filename = 'output/%s/cluster_properties.csv' % cluster
    reader = {}
    if os.path.isfile(filename):
        with open(filename, 'r') as F:
            reader = csv.DictReader(F)
            reader = [dict(row) for row in reader][0]

    reader[key] = value

    with open(filename, 'w') as F:
        writer = csv.DictWriter(F, reader.keys())
        writer.writeheader()
        writer.writerow(reader)
