from .observations import ListDates
from .astrometry import GetAstrometryOffset
from .spectral_type import AbsMag, SpectralType
from .gaia import GetRParams, GetGaiaOutliers
from .io import Binning, GetWCS, GetFITSValues, WriteClusterProperty
from .be_filter import GetVMagLimit

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.gridspec import GridSpec

import numpy as np
from astropy.io import fits

import os.path
import csv
from collections import OrderedDict


class Analysis:

    def __init__(self, cluster, app):
        self.cluster = cluster
        self.app = app

        self.BeCandidates = []
        self.rej_BeCandidates = []
        self.outlier_BeCandidates = []

        self.mostB_date = None

        self.distance_mean, self.distance_std = GetRParams(self.cluster)

    def ProcessAnalysis(self):
        """Controls full analysis process.

        This program compiles a full list of Be candidates from the calculated Be
        candidates of individual nights, finds corresponding targets from other
        dates, then summarizes the results in output.

        Args:
            cluster (str): Cluster for which analysis will be conducted.
            app (Application): The GUI application object that controls processing.

        """
        # Compile list of Be stars
        self.CompileBeLists()
        for ls in ('belist', 'rej_belist'):
            self.FindCorrespondingTargets(ls)

        # Output
        self.NightSummary()

        for ls in ('belist', 'rej_belist'):
            GetExcess(self, ls)

        self.SortBelist()
        for ls in ('belist', 'rej_belist', 'outlier_belist'):
            self.BeSummary(ls)

        self.ColorAnalysis()
        self.SpectralTypeDist()

        for date in ListDates(self.cluster):
            self.BeCandidatePlots(date)

    class BeCandidate:

        def __init__(self):
            self.x = None
            self.y = None
            self.xref = None
            self.yref = None
            self.ra = None
            self.dec = None
            self.julian = None
            self.date = None
            self.count = None

            self.phot = {}
            self.excess = None
            self.excess_err = None

            self.absmag = None
            self.spectraltype = None

            self.observed_be = True

        def GetAbsMag(self, distance):
            if self.absmag is None:
                self.absmag = AbsMag(self.phot['Vmag'], distance)

            return self.absmag

        def GetSpectralType(self, distance):
            if self.spectraltype is None:
                self.spectraltype = SpectralType(self.phot['Vmag'], distance)

            return self.spectraltype

    def MatchTarget_CBL(self, tol, coords, data):
        potentialMatches = []
        for target in data:
            if abs(coords[0] - target.xref) <= tol and \
               abs(coords[1] - target.yref) <= tol:
                potentialMatches.append(target)

        # Match with the closest target
        if potentialMatches:
            d = [np.sqrt((coords[0] - x.xref)**2 + (coords[1] - x.yref)**2)
                 for x in potentialMatches]
            d = np.array(d)
            match = potentialMatches[np.argmin(d)]
            return match

        return None

    def CompileBeLists(self):
        """Creates data sets for all Be candidates from each date of observation.

        This reads from output from each night to compile each indivual Be candidate
        and uniquely identify them.  Targets are then put into an 'accepted' or
        'rejected' list based on calculated cluster membership.

        Args:
            cluster (str): Cluster for which Be candidates will be compiled.
            app (Application): The GUI application object that controls processing.

        """
        acc_count = 1
        rej_count = 1
        baseDate = ListDates(self.cluster)[0]

        for date in ListDates(self.cluster):
            print("  Retrieving Be data from %s..." % date)

            filename = 'output/%s/%s/belist_scaled.dat' % (self.cluster, date)
            data = np.loadtxt(filename, ndmin=2)

            # Get header info
            fits = GetFITSValues(self.cluster, date,
                                 ['JD', 'XBINNING', 'RA', 'DEC'])
            julian = fits['JD']
            binning = fits['XBINNING']

            # Get ra/decs of CERTAIN outliers (by distance) for this date
            outliers = GetGaiaOutliers(self.cluster, date)
            w = GetWCS(self.cluster, date)

            # Get x- and y-offsets
            if date != baseDate:
                xOffset, yOffset = GetAstrometryOffset(self.cluster, date,
                                                       baseDate)
            else:
                xOffset, yOffset = (0, 0)

            for target in data:
                # Compare with other dates, see if target is already found
                xref = binning * target[0] + xOffset
                yref = binning * target[1] + yOffset
                match = self.MatchTarget_CBL(self.app.cooTol, (xref, yref),
                                             self.BeCandidates + self.rej_BeCandidates)
                if match is not None:
                    # If target refers to an already listed star
                    ra = match.ra
                    dec = match.dec
                    count = match.count
                    accept = match in self.BeCandidates
                else:
                    # If target does not match with a listed star
                    radec = w.all_pix2world(target[0], target[1], 0)
                    ra, dec = tuple([float(i) for i in radec])

                    if (ra, dec) not in outliers:
                        count = acc_count
                        acc_count += 1
                        accept = True
                    else:
                        count = rej_count
                        rej_count += 1
                        accept = False

                be = self.BeCandidate()

                be.x = target[0]
                be.y = target[1]
                be.xref = xref
                be.yref = yref
                be.ra = ra
                be.dec = dec
                be.julian = julian
                be.date = date
                be.count = count

                phots = ('Bmag', 'Berr_low', 'Berr_high',
                         'Vmag', 'Verr_low', 'Verr_high',
                         'Rmag', 'Rerr_low', 'Rerr_high',
                         'Hmag', 'Herr_low', 'Herr_high')
                for p, i in zip(phots, range(2, 14)):
                    be.phot[p] = target[i]

                be.observed_be = True

                if accept:
                    self.BeCandidates.append(be)
                else:
                    self.rej_BeCandidates.append(be)

    def MatchTarget_FCT(self, tol, coords, phot):
        potentialMatches = []

        for target in phot:
            if abs(coords[0] - target[14]) <= tol and \
               abs(coords[1] - target[15]) <= tol:
                potentialMatches.append(target)

        # Match with the closest target
        if potentialMatches:
            d = [np.sqrt((coords[0] - x[14])**2 + (coords[1] - x[15])**2)
                 for x in potentialMatches]
            d = np.array(d)
            match = potentialMatches[np.argmin(d)]
            return match

        return None

    def FindCorrespondingTargets(self, belist):
        """Finds targets from other nights that correspond to Be candidates.

        Uses matching to determine which targets from each night of observation
        identify with each star in given Be candidate list.

        Args:
            cluster (str): Cluster for which Be candidates will be compiled.
            app (Application): The GUI application object that controls processing.
            belist (list): List (accepted or rejected) to check for corresponding
                           targets.

        """
        dates = ListDates(self.cluster)
        baseDate = dates[0]

        if belist == 'belist':
            belist = self.BeCandidates
        else:
            belist = self.rej_BeCandidates

        max_count = max((x.count for x in belist))

        for date in dates:
            filename = 'output/%s/%s/phot_scaled.dat' % (self.cluster, date)
            phot = np.loadtxt(filename).tolist()

            # Get header info
            fits = GetFITSValues(self.cluster, date, ['JD', 'XBINNING'])
            julian = fits['JD']
            binning = fits['XBINNING']

            # Get x- and y-offsets
            if date != baseDate:
                xOffset, yOffset = GetAstrometryOffset(self.cluster, date,
                                                       baseDate)
            else:
                xOffset, yOffset = (0, 0)

            for count in range(1, max_count + 1):
                candidate_full = [x for x in belist if x.count == count]

                if date in [x.date for x in candidate_full]:
                    continue

                xref = np.mean([x.xref for x in candidate_full])
                yref = np.mean([x.yref for x in candidate_full])

                ra = candidate_full[0].ra
                dec = candidate_full[0].dec

                for target in phot:
                    target.append(binning * target[0] + xOffset)
                    target.append(binning * target[1] + yOffset)

                match = self.MatchTarget_FCT(self.app.cooTol, (xref, yref), phot)

                if match is None:
                    continue

                phot.remove(match)

                be = self.BeCandidate()

                be.x = match[0]
                be.y = match[1]
                be.xref = match[14]
                be.yref = match[15]
                be.ra = ra
                be.dec = dec
                be.julian = julian
                be.date = date
                be.count = count

                phots = ('Bmag', 'Berr_low', 'Berr_high',
                         'Vmag', 'Verr_low', 'Verr_high',
                         'Rmag', 'Rerr_low', 'Rerr_high',
                         'Hmag', 'Herr_low', 'Herr_high')
                for p, i in zip(phots, range(2, 14)):
                    be.phot[p] = match[i]

                be.observed_be = False

                belist.append(be)

    def GetBeRatio(self, be, b):
        try:
            ratio = be / (be + b)
            error = np.sqrt(be / (b + be) * (1 - be / (b + be)) / (b + be))
        except ZeroDivisionError:
            ratio = 0
            error = 0

        return (ratio, error)

    def FindBStars(self, date, rejected=False):
        """Finds B-type stars without H-alpha excess for given date.

        Using photometric bounds and calculated Be candidates as well as cluster
        membership calculations, B-type targets are extracted.

        Args:
            cluster (str): Cluster for which B stars will be found.
            app (Application): The GUI application object that controls processing.
            date (str): Date for which B stars will be found.
            BeCandidates (list): List of accepted Be candidates.
            rej_BeCandidates (list): List of rejected Be candidates.
            rejected (bool): Determines if accepted or rejected B-type stars are
                             desired. False if accepted is desired, and vice versa.

        Returns:
            list: List of B-type targets

        """
        if not rejected:
            filename = 'output/%s/%s/phot_scaled_accepted.dat' % (self.cluster, date)
        else:
            filename = 'output/%s/%s/phot_scaled_rejected.dat' % (self.cluster, date)

        data = np.loadtxt(filename).tolist()
        BData = []

        # Find all B- and Be-type stars
        Vlim = GetVMagLimit(self.cluster)
        for target in data:
            b_v = target[2] - target[5]
            if self.app.B_VMin < b_v < self.app.B_VMax and target[5] < Vlim:
                BData.append(target)

        # Pick out observed Be stars to result in ONLY B-type stars
        baseDate = ListDates(self.cluster)[0]
        xOffset, yOffset = GetAstrometryOffset(self.cluster, date, baseDate)
        binning = Binning(self.cluster, date)

        for b in reversed(BData):
            xref = b[0] * binning + xOffset
            yref = b[1] * binning + yOffset
            for be in self.BeCandidates + self.rej_BeCandidates:
                if abs(xref - be.xref) <= self.app.cooTol and \
                   abs(yref - be.yref) <= self.app.cooTol:
                    BData.remove(b)
                    break

        return BData

    def NightSummary(self):
        """Writes a summary of each night to a file.

        Returns:
            str: Returns the date with the most B-type targets observed.

        """
        filename = 'output/%s/summary_dates.txt' % self.cluster
        summary_file = open(filename, 'w')

        t = '================================================\n' + \
            '                %s Summary                  \n' % self.cluster + \
            '================================================\n\n\n'
        summary_file.write(t)

        NumBDict = {}

        for date in ListDates(self.cluster):
            t = '                   %s                   \n' % date + \
                '------------------------------------------------\n'
            summary_file.write(t)

            # Technical information
            filename = 'photometry/%s/%s/B1.fits' % (self.cluster, date)
            binning = Binning(self.cluster, date)
            summary_file.write('Binning: %d\n' % binning)

            summary_file.write('Exposures: ')
            files = ['B1', 'B3', 'V1', 'V3', 'R1', 'R3', 'H1', 'H3']
            for file in files:
                filename = 'photometry/%s/%s/%s.fits' % (self.cluster, date, file)
                G = fits.getheader(filename)
                summary_file.write('%s: %d' % (file, G['EXPTIME']))
                if file != files[len(files) - 1]:
                    summary_file.write(' | ')

            summary_file.write('\n\n')

            # Threshold information
            filename = 'output/%s/%s/thresholds.dat' % (self.cluster, date)
            thresholds = np.loadtxt(filename)
            t = "Thresholds:\n" + \
                "   Constant: %.3f     Linear: %.3f, %.3f\n\n" % \
                (thresholds[0][1], thresholds[1][0], thresholds[1][1])

            summary_file.write(t)

            # Be candidate and B star information
            filename = 'output/%s/%s/belist_scaled.dat' % (self.cluster, date)
            beData = np.loadtxt(filename)
            NumBe = len(beData)

            bData = self.FindBStars(date)
            NumB = len(bData)

            be_ratio = self.GetBeRatio(NumBe, NumB)

            t = "Be candidates:              %d\n" % NumBe + \
                "B stars:                    %d\n" % NumB + \
                "Be ratio:                   %.3f +/- %.3f\n\n\n" % be_ratio

            summary_file.write(t)

            # Individual date info
            filename = 'output/%s/%s/count_and_ratio.txt' % (self.cluster, date)
            with open(filename, 'w') as F:
                F.write('%d\n' % NumBe)
                F.write('%.3f\n' % be_ratio[0])
                F.write('%.3f' % be_ratio[1])

            NumBDict[date] = NumB

        summary_file.close()
        self.mostB_date = max(NumBDict, key=NumBDict.get)

    def GetExcess(self, belist):
        if belist == 'belist':
            data = self.BeCandidates
        elif belist == 'rej_belist':
            data = self.rej_BeCandidates

        for target in data:
            filename = 'output/%s/%s/thresholds.dat' % (self.cluster, target.date)
            thresholds = np.loadtxt(filename)
            if self.app.threshold_type == 'Constant':
                slope = thresholds[0][0]
                intercept = thresholds[0][1]
            elif self.app.threshold_type == 'Linear':
                slope = thresholds[1][0]
                intercept = thresholds[1][1]

            r_h = target.phot['Rmag'] - target.phot['Hmag']
            b_v = target.phot['Bmag'] - target.phot['Vmag']
            r_h0 = slope * b_v + intercept
            target.excess = r_h - r_h0

            target.excess_err = np.sqrt(target.phot['Rerr_low']**2 +
                                        target.phot['Herr_high']**2)

    def SeparateOutliers(self):
        max_count = max((x.count for x in self.BeCandidates))
        for i in range(1, max_count + 1):
            candidate_full = [x for x in self.BeCandidates if x.count == i]
            if len([x for x in candidate_full if x.excess > 0]) == 1:
                for target in candidate_full:
                    self.BeCandidates.remove(target)
                    self.outlier_BeCandidates.append(target)

    def SortBelist(self):
        # Sort by identifier (count) then date
        self.BeCandidates = sorted(self.BeCandidates,
                                   key=lambda x: (x.count, x.date))
        self.rej_BeCandidates = sorted(self.rej_BeCandidates,
                                       key=lambda x: (x.count, x.date))
        self.outlier_BeCandidates = sorted(self.outlier_BeCandidates,
                                           key=lambda x: (x.count, x.date))

    def BeSummary(self, belist):
        """Creates finalized output for Be candidate lists.

        This finalized summary includes the determination of each target's spectral
        type and its excess R-H value from the calculated threshold for each night.
        Also analyzes spectral type distribution through histograms and ratios.

        Args:
            belist (str): Specific list of Be candidates to analyze.

        """
        if belist == 'belist':
            data = self.BeCandidates
            filename = 'output/%s/BeList.txt' % self.cluster
        elif belist == 'rej_belist':
            data = self.rej_BeCandidates
            filename = 'output/%s/BeList_rejected.txt' % self.cluster
        elif belist == 'outlier_belist':
            data = self.outlier_BeCandidates
            filename = 'output/%s/BeList_outliers.txt' % self.cluster

        if not data:
            return

        with open(filename, 'w') as file:
            for target in data:
                # Write to file
                t = '%s-WBBe%d' % (self.cluster, target.count) + '\t' + \
                    '%.3f' % target.GetAbsMag(self.distance_mean) + '\t' + \
                    target.GetSpectralType(self.distance_mean) + '\t' + \
                    '%.10f' % target.ra + '\t' + \
                    '%.10f' % target.dec + '\t' + \
                    '%.10f' % target.julian + '\t' + \
                    '%.3f' % target.phot['Bmag'] + '\t' + \
                    '%.3f' % target.phot['Berr_low'] + '\t' + \
                    '%.3f' % target.phot['Berr_high'] + '\t' + \
                    '%.3f' % target.phot['Vmag'] + '\t' + \
                    '%.3f' % target.phot['Verr_low'] + '\t' + \
                    '%.3f' % target.phot['Verr_high'] + '\t' + \
                    '%.3f' % target.phot['Rmag'] + '\t' + \
                    '%.3f' % target.phot['Rerr_low'] + '\t' + \
                    '%.3f' % target.phot['Rerr_high'] + '\t' + \
                    '%.3f' % target.phot['Hmag'] + '\t' + \
                    '%.3f' % target.phot['Herr_low'] + '\t' + \
                    '%.3f' % target.phot['Herr_high'] + '\t' + \
                    '%.3f' % target.excess + '\t' + \
                    '%.3f' % target.excess_err + '\n'

                file.write(t)

    def ColorAnalysis(self):
        path = 'output/' + self.cluster + '/color_color_plots'
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            for file in os.listdir(path):
                file_path = os.path.join(path, file)
                if os.path.isfile(file_path):
                    os.remove(file_path)

        filename = 'output/%s/BeList_colors.txt' % self.cluster
        file = open(filename, 'w')

        plt.style.use('researchpaper')
        fig = plt.figure(dpi=200)
        fig.set_size_inches(16, 11)

        num_cols = 4
        gs = GridSpec(3, num_cols)

        max_count = max((x.count for x in self.BeCandidates))
        page = 1

        date_to_color_map = {
            '201505': ('#ea8f34', '#ffbf87'),   # orange
            '201508': ('#ea3535', '#ff8787'),   # red
            '201511': ('#b933ea', '#c877e5'),   # purple
            '201512': ('#028ecc', '#5ec3f2'),   # blue
            '201610': ('#13a019', '#57d65d')    # green
        }

        for i in range(1, max_count + 1):
            candidate_full = [x for x in self.BeCandidates if x.count == i]
            identifier = '%s-WBBe%d' % (self.cluster, i)

            col = (i - 1) % num_cols
            ax_vr = fig.add_subplot(gs[0, col])
            ax_br = fig.add_subplot(gs[1, col])
            ax_rh = fig.add_subplot(gs[2, col])

            b_v = []
            b_v_err_low = []
            b_v_err_high = []
            v_r = []
            v_r_err_low = []
            v_r_err_high = []
            b_r = []
            b_r_err_low = []
            b_r_err_high = []
            r_h = []
            r_h_err_low = []
            r_h_err_high = []
            excess_shown = []
            dates = []

            for target in candidate_full:
                bv = target.phot['Bmag'] - target.phot['Vmag']
                vr = target.phot['Vmag'] - target.phot['Rmag']
                br = target.phot['Bmag'] - target.phot['Rmag']
                rh = target.phot['Rmag'] - target.phot['Hmag']

                b_v.append(bv)
                v_r.append(vr)
                b_r.append(br)
                r_h.append(rh)

                bv_err_low = np.sqrt(target.phot['Berr_low']**2 +
                                     target.phot['Verr_high']**2)
                vr_err_low = np.sqrt(target.phot['Verr_low']**2 +
                                     target.phot['Rerr_high']**2)
                br_err_low = np.sqrt(target.phot['Berr_low']**2 +
                                     target.phot['Rerr_high']**2)
                rh_err_low = np.sqrt(target.phot['Rerr_low']**2 +
                                     target.phot['Herr_high']**2)

                b_v_err_low.append(bv_err_low)
                v_r_err_low.append(vr_err_low)
                b_r_err_low.append(br_err_low)
                r_h_err_low.append(rh_err_low)

                bv_err_high = np.sqrt(target.phot['Berr_high']**2 +
                                      target.phot['Verr_low']**2)
                vr_err_high = np.sqrt(target.phot['Verr_high']**2 +
                                      target.phot['Rerr_low']**2)
                br_err_high = np.sqrt(target.phot['Berr_high']**2 +
                                      target.phot['Rerr_low']**2)
                rh_err_high = np.sqrt(target.phot['Rerr_high']**2 +
                                      target.phot['Herr_low']**2)

                b_v_err_high.append(bv_err_high)
                v_r_err_high.append(vr_err_high)
                b_r_err_high.append(br_err_high)
                r_h_err_high.append(rh_err_high)

                dates.append(target.date[:6])
                excess_shown.append(target.observed_be)

                # Write to file
                t = identifier + '\t' + \
                    '%.10f' % target.julian + '\t' + \
                    '%.3f' % bv + '\t' + \
                    '%.3f' % vr + '\t' + \
                    '%.3f' % br + '\t' + \
                    '%.3f' % rh + '\n'

                file.write(t)

            # Plot lines
            for j in range(len(b_v) - 1):
                if dates[j] != dates[j + 1]:
                    continue

                line_bv = (b_v[j], b_v[j + 1])
                line_vr = (v_r[j], v_r[j + 1])
                line_br = (b_r[j], b_r[j + 1])
                line_rh = (r_h[j], r_h[j + 1])

                ax_vr.plot(line_bv, line_vr, '-', color='#3a3a3a', zorder=1)
                ax_br.plot(line_bv, line_br, '-', color='#3a3a3a', zorder=1)
                ax_rh.plot(line_bv, line_rh, '-', color='#3a3a3a', zorder=1)

                mean_bv = np.mean(line_bv)
                ax_vr.annotate('', xytext=(line_bv[0], line_vr[0]),
                               xy=(mean_bv, np.mean(line_vr)),
                               arrowprops=dict(arrowstyle="-|>", color='#3a3a3a'),
                               size=18, zorder=1)
                ax_br.annotate('', xytext=(line_bv[0], line_br[0]),
                               xy=(mean_bv, np.mean(line_br)),
                               arrowprops=dict(arrowstyle="-|>", color='#3a3a3a'),
                               size=18, zorder=1)
                ax_rh.annotate('', xytext=(line_bv[0], line_rh[0]),
                               xy=(mean_bv, np.mean(line_rh)),
                               arrowprops=dict(arrowstyle="-|>", color='#3a3a3a'),
                               size=18, zorder=1)

            # Plot markers
            s = 100
            for bv, vr, br, rh, \
                bv_err_low, bv_err_high, vr_err_low, vr_err_high, \
                br_err_low, br_err_high, rh_err_low, rh_err_high, \
                ex, date in \
                    zip(b_v, v_r, b_r, r_h,
                        b_v_err_low, b_v_err_high, v_r_err_low, v_r_err_high,
                        b_r_err_low, b_r_err_high, r_h_err_low, r_h_err_high,
                        excess_shown, dates):
                color = date_to_color_map[date]
                if ex:
                    mk = '^'
                    facecolor = color[0]
                    lw = None
                else:
                    mk = 'o'
                    facecolor = 'none'
                    lw = 2
                ax_vr.scatter(bv, vr, marker=mk,
                              facecolors=facecolor, edgecolors=color[0],
                              linewidths=lw, s=s, zorder=3)
                ax_br.scatter(bv, br, marker=mk,
                              facecolors=facecolor, edgecolors=color[0],
                              linewidths=lw, s=s, zorder=3)
                ax_rh.scatter(bv, rh, marker=mk,
                              facecolors=facecolor, edgecolors=color[0],
                              linewidths=lw, s=s, zorder=3)

                ax_vr.errorbar(bv, vr, xerr=[[bv_err_low], [bv_err_high]],
                               yerr=[[vr_err_low], [vr_err_high]], fmt='none',
                               ecolor=color[1], elinewidth=2, zorder=2)
                ax_br.errorbar(bv, br, xerr=[[bv_err_low], [bv_err_high]],
                               yerr=[[br_err_low], [br_err_high]], fmt='none',
                               ecolor=color[1], elinewidth=2, zorder=2)
                ax_rh.errorbar(bv, rh, xerr=[[bv_err_low], [bv_err_high]],
                               yerr=[[rh_err_low], [rh_err_high]], fmt='none',
                               ecolor=color[1], elinewidth=2, zorder=2)

            ax_vr.set_title(identifier, fontsize=24)

            fontsize = 20
            ax_rh.set_xlabel('B-V', fontsize=fontsize)

            if col == 0:
                ax_vr.set_ylabel('V-R', fontsize=fontsize)
                ax_br.set_ylabel('B-R', fontsize=fontsize)
                ax_rh.set_ylabel('R-H', fontsize=fontsize)

            # label_s = '%.2f'
            # label_maj_x = 0.02
            for ax in (ax_vr, ax_br, ax_rh):
                ax.tick_params(axis='both', which='major', labelsize=14, width=1,
                               length=5)
                # ax.xaxis.set_major_locator(MultipleLocator(label_maj_x))
                # ax.xaxis.set_major_formatter(FormatStrFormatter(label_s))
                # ax.yaxis.set_major_formatter(FormatStrFormatter(label_s))

            if i % num_cols != 0 and i != max_count:
                continue

            # Figure legend
            dates_to_formal = {
                '201505': 'May. 2015',
                '201508': 'Aug. 2015',
                '201511': 'Nov. 2015',
                '201512': 'Dec. 2015',
                '201610': 'Oct. 2016'
            }

            fc = 'none'
            if excess_shown[0]:
                mk = '^'
                lw = None
            else:
                mk = 'o'
                lw = 2

            dates_unique = list(OrderedDict.fromkeys(dates))
            markers = []
            for date in dates_unique:
                color = date_to_color_map[date][0]
                if excess_shown[0]:
                    fc = color
                d, = ax_vr.plot(b_v[0], v_r[0], '^', markeredgecolor=color,
                                markeredgewidth=lw, markerfacecolor=fc)
                markers.append(d)

            dates_unique = [dates_to_formal[x] for x in dates_unique]
            fig.legend(markers, dates_unique, loc='right',
                       bbox_to_anchor=(1.15, 0.5), markerscale=1.2,
                       fontsize=18, frameon=False)

            # Output
            plt.tight_layout()

            filename = 'output/%s/color_color_plots/2cd_%d.png' % (self.cluster, page)
            fig.savefig(filename, dpi=200, bbox_inches="tight")
            page += 1
            plt.clf()

        file.close()

    def SpectralTypeDist(self):
        BData = self.FindBStars(self.mostB_date)

        # Be ratios by spectral type
        be_spectralTypes = []
        count = (x.count for x in self.BeCandidates)
        for i in range(1, max(count) + 1):
            m = [x.phot['Vmag'] for x in self.BeCandidates if x.count == i]
            s = SpectralType(np.mean(m), self.distance_mean)
            be_spectralTypes.append(s)

        b_spectralTypes = [SpectralType(target[5], self.distance_mean) for target in BData]

        frequencies = []
        spec_types = {
            'Unknown': ('--'),
            'O6-O8': ('O6', 'O7', 'O8'),
            'B0-B3': ('B0', 'B1', 'B2', 'B3'),
            'B4-B5': ('B4', 'B5'),
            'B6-B9': ('B6', 'B7', 'B8', 'B9'),
            'A': ('A0', 'A1', 'A2', 'A5', 'A8')
        }

        filename = 'output/%s/ratios_spec_type.csv' % self.cluster
        with open(filename, 'w') as F:
            writer = csv.writer(F)
            header = ('spec_type', 'no_Be', 'no_B', 'ratio', 'ratio_err')
            writer.writerow(header)

            for spec_type in spec_types:
                no_be = len([x for x in be_spectralTypes if x in spec_types[spec_type]])
                no_b = len([x for x in b_spectralTypes if x in spec_types[spec_type]])
                be_ratio, be_ratio_err = self.GetBeRatio(no_be, no_b)
                writer.writerow((spec_type, '%d' % no_be, '%d' % no_b,
                                 '%.3f' % be_ratio, '%.3f' % be_ratio_err))
                frequencies.append(no_be)

        # Plot spectral type histogram
        plt.style.use('researchpaper')
        fig, ax = plt.subplots()

        x_coordinates = np.arange(len(frequencies))
        ax.bar(x_coordinates, frequencies, align='center', color='#3f3f3f')
        ax.set_xlabel('Spectral Type')
        ax.set_ylabel('Frequency')

        ax.xaxis.set_major_locator(plt.FixedLocator(x_coordinates))
        ax.xaxis.set_major_formatter(plt.FixedFormatter(list(spec_types.keys())))

        spine_lw = 4
        [ax.spines[axis].set_linewidth(spine_lw)
         for axis in ['top', 'bottom', 'left', 'right']]

        filename = 'output/%s/spectral_type_dist.png' % self.cluster
        fig.savefig(filename)

        plt.close('all')

        # Be ratio and percentages for full cluster
        no_be = max(x.count for x in self.BeCandidates)
        no_b = len(BData)
        be_ratio = self.GetBeRatio(no_be, no_b)

        stables = []

        for k in range(1, no_be + 1):
            candidates = [x for x in self.BeCandidates if x.count == k]
            missed_lim = 0
            times_missed = len([x for x in candidates if not x.observed_be])

            if times_missed <= missed_lim:
                stables.append(k)

        no_stables = len(stables)
        no_transients = no_be - no_stables

        percentage = no_transients / no_be
        percentage_err = np.sqrt(percentage * (1 - percentage) / no_be)

        WriteClusterProperty(self.cluster, 'no_Be', '%d' % no_be)
        WriteClusterProperty(self.cluster, 'ratio', '%.3f' % be_ratio[0])
        WriteClusterProperty(self.cluster, 'ratio_err', '%.3f' % be_ratio[1])
        WriteClusterProperty(self.cluster, 'pct', '%.3f' % percentage)
        WriteClusterProperty(self.cluster, 'pct_err', '%.3f' % percentage_err)

    def BeCandidatePlots(self, date):
        """Creates detailed CMDs and 2CDs for a specified date.

        Detailed plots created are CMDs and 2CDs that specify cluster membership
        of all observed Be candidates and B-type stars.

        Args:
            cluster (str): Cluster for which information is written.
            app (Application): The GUI application object that controls processing.
            date (str): Date for which information is written.
            BeCandidates (list): List of accepted Be candidates.
            rej_BeCandidates (list): List of rejected Be candidates.

        """
        filename = 'output/%s/%s/phot_scaled.dat' % (self.cluster, date)
        phot = np.loadtxt(filename, ndmin=2).tolist()

        # Be candidates in cluster
        Be_in = []
        for candidate in (x for x in self.BeCandidates if x.date == date):
            for target in phot:
                if target[0] == candidate.x and target[1] == candidate.y:
                    Be_in.append(target)
                    phot.remove(target)
                    break

        # Be candidates outside cluster
        Be_out = []
        for candidate in (x for x in self.rej_BeCandidates if x.date == date):
            for target in phot:
                if target[0] == candidate.x and target[1] == candidate.y:
                    Be_out.append(target)
                    phot.remove(target)
                    break

        # B-type stars in cluster
        B_in = self.FindBStars(date)

        # B-type stars outside cluster
        B_out = self.FindBStars(date, rejected=True)

        # Plot data
        filename = 'output/%s/%s/phot_scaled.dat' % (self.cluster, date)
        data = np.loadtxt(filename)

        B_V = data[:, 2] - data[:, 5]
        B_Verr_low = np.sqrt(data[:, 3]**2 + data[:, 7]**2)
        B_Verr_high = np.sqrt(data[:, 4]**2 + data[:, 6]**2)
        R_H = data[:, 8] - data[:, 11]
        R_Herr_low = np.sqrt(data[:, 9]**2 + data[:, 13]**2)
        R_Herr_high = np.sqrt(data[:, 10]**2 + data[:, 12]**2)

        def Plot(plot_type):
            plt.style.use('researchpaper')
            fig, ax = plt.subplots()

            if plot_type == 'cmd':
                # Plotting
                ax.plot(B_V, data[:, 5], 'o', color='#3d3d3d', markersize=10)
                ax.errorbar(B_V, data[:, 5], xerr=[B_Verr_low, B_Verr_high],
                            yerr=[data[:, 6], data[:, 7]], fmt='none',
                            ecolor='#8c8c8c', elinewidth=7)

                ax.plot([x[2] - x[5] for x in B_out], [x[5] for x in B_out],
                        'o', color='#6ba3ff', markersize=12, markeredgewidth=2,
                        label='B outside cluster')
                ax.plot([x[2] - x[5] for x in B_in], [x[5] for x in B_in],
                        'o', color='#ff5151', markersize=12, markeredgewidth=2,
                        label='B in cluster')
                ax.plot([x[2] - x[5] for x in Be_out], [x[5] for x in Be_out],
                        'D', color='#6ba3ff', markersize=12, markeredgewidth=2,
                        markeredgecolor='#2d2d2d', label='Be outside cluster')
                ax.plot([x[2] - x[5] for x in Be_in], [x[5] for x in Be_in],
                        'D', color='#ff5151', markersize=12, markeredgewidth=2,
                        markeredgecolor='#2d2d2d', label='Be in cluster')

                # Settings
                ax.set_ylabel('V')
                ax.invert_yaxis()

                ax.yaxis.set_major_locator(MultipleLocator(2))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
                ax.yaxis.set_minor_locator(MultipleLocator(0.5))

                output = 'CMD_detailed'

            elif plot_type == '2cd':
                # Plotting
                ax.plot(B_V, R_H, 'o', color='#3d3d3d', markersize=10)
                ax.errorbar(B_V, R_H, xerr=[B_Verr_low, B_Verr_high],
                            yerr=[R_Herr_low, R_Herr_high], fmt='none',
                            ecolor='#8c8c8c', elinewidth=7)

                ax.plot([x[2] - x[5] for x in B_out],
                        [x[8] - x[11] for x in B_out], 'o', color='#6ba3ff',
                        markersize=12, markeredgewidth=2,
                        label='B outside cluster')
                ax.plot([x[2] - x[5] for x in B_in],
                        [x[8] - x[11] for x in B_in], 'o', color='#ff5151',
                        markersize=12, markeredgewidth=2,
                        label='B in cluster')
                ax.plot([x[2] - x[5] for x in Be_out],
                        [x[8] - x[11] for x in Be_out], 'D', color='#6ba3ff',
                        markersize=12, markeredgewidth=2,
                        markeredgecolor='#2d2d2d', label='Be outside cluster')
                ax.plot([x[2] - x[5] for x in Be_in],
                        [x[8] - x[11] for x in Be_in], 'D', color='#ff5151',
                        markersize=12, markeredgewidth=2,
                        markeredgecolor='#2d2d2d', label='Be in cluster')

                # Settings
                ax.set_ylabel('R-Halpha')

                ax.yaxis.set_major_locator(MultipleLocator(1))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
                ax.yaxis.set_minor_locator(MultipleLocator(0.25))

                filename = 'output/%s/%s/thresholds.dat' % (self.cluster, date)
                thresholds = np.loadtxt(filename)
                if self.app.threshold_type == 'Constant':
                    file = thresholds[0]
                elif self.app.threshold_type == 'Linear':
                    file = thresholds[1]
                slope = file[0]
                intercept = file[1]

                linex = np.array([self.app.B_VMin, self.app.B_VMax])
                liney = slope * linex + intercept
                ax.plot(linex, liney, '--', color='#ff5151',
                        label='Be Threshold', linewidth=6)

                output = '2CD_detailed'

            # Shared settings
            ax.set_xlabel('B-V')

            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax.xaxis.set_minor_locator(MultipleLocator(0.25))

            spine_lw = 4
            [ax.spines[axis].set_linewidth(spine_lw)
             for axis in ['top', 'bottom', 'left', 'right']]

            ax.legend(fontsize=20)

            # Output
            filename = 'output/%s/%s/plots/%s.png' % (self.cluster, date, output)
            fig.savefig(filename)

            # fig.patch.set_facecolor('#f4f4f4')
            ax.set_facecolor('#f4f4f4')
            filename = 'output/%s/%s/plots/%s_web.png' % (self.cluster, date, output)
            fig.savefig(filename, facecolor='#f4f4f4')

            plt.close('all')

        Plot('2cd')
        Plot('cmd')
