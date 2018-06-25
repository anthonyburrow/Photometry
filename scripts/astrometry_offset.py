import os
from numpy import mean
import simplejson
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.application import MIMEApplication
from email.encoders import encode_noop
from my_generator import MyGenerator
from io import StringIO
import time
from astropy.io import fits
from observations import Observations


class AstrometryOffset:
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
        baseDate = Observations.ListDates(cluster)[0]

        # Create reference plate scaled image
        base_wcs = "../photometry/" + cluster + "/" + baseDate + "/B1_wcs.fits"
        if not os.path.isfile(base_wcs):
            upload = "../photometry/" + cluster + "/" + baseDate + "/B1.fits"
            self.GetPlateScaledFile(upload, base_wcs)

        # Create plate scaled image for specified date
        upload = "photometry/" + cluster + "/" + date + "/"
        self.GetPlateScaledFile(upload + "B1.fits", upload + "B1_wcs.fits")

        # Get accurate RA and DEC values
        # with fits.open("../photometry/" + cluster + "/" + baseDate + "/") as file:
        #     baseRA = file[0].header["CRVAL1"]
        #     baseDec = file[0].header["CRVAL2"]
        #
        # with fits.open("../photometry/" + cluster + "/" + date + "/") as file:
        #     RA = file[0].header["CRVAL1"]
        #     Dec = file[0].header["CRVAL2"]

        # Get new information from astrometry.net
        with fits.open("../photometry/" + cluster + "/" + baseDate + "/B1.fits") as file:
            baseMaxPixels = file[0].header['NAXIS1']
            baseBin = file[0].header['XBINNING']
        with fits.open("../photometry/" + cluster + "/" + date + "/B1.fits") as file:
            Bin = file[0].header['XBINNING']

        with fits.open("../photometry/" + cluster + "/" + baseDate + "/B1_corr.fits") as file:
            baseCorr = file[1].data
        with fits.open("../photometry/" + cluster + "/" + date + "/B1_corr.fits") as file:
            corr = file[1].data

        # Calculate offsets
        count = 0
        sample = []
        for baseTarget in baseCorr:
            if baseTarget[0] >= baseMaxPixels * 0.25 and baseTarget[0] <= baseMaxPixels * 0.75 and \
               baseTarget[1] >= baseMaxPixels * 0.25 and baseTarget[1] <= baseMaxPixels * 0.75:
                for target in corr:
                    if baseTarget[6] == target[6] and baseTarget[7] == target[7]:
                        xOff = baseBin * baseTarget[0] - Bin * target[0]
                        yOff = baseBin * baseTarget[1] - Bin * target[1]
                        sample.append([xOff, yOff])
            if count == 5:
                break

        offsets = [mean(sample[:, 0]), mean(sample[:, 1])]
        return offsets

    def GetPlateScaledFile(self, upload_filename, output_path):
        """Get the plate scaled information from astrometry.net

        Uploads a specified .fits file and retrieves output information from astrometry.net.  Uses the
        technique provided by API documetation and the contribution provided by Christopher Klein.

        Args:
            upload_filename: The filename (and path) of the input .fits image.
            ouput_path: The directory that the output is to be placed.

        """

        # Get API key
        if not os.path.exists("etc/"):
            os.makedirs("etc/")
        key = "etc/apikey"
        if os.path.isfile(key):
            with open(key) as F:
                apikey = F.readline()
        else:
            with open(key, 'w') as F:
                F.write(input("Astrometry.net API key (found on profile): "))

        # Login to the service
        apiurl = "http://nova.astrometry.net/api/"
        login_url = apiurl + "login"
        login_args = {'apikey': apikey}
        json = simplejson.dumps(login_args)
        login_data = {'request-json': json}
        login_data = urlencode(login_data)
        login_headers = {}
        request = Request(url=login_url, headers=login_headers, data=login_data)

        f = urlopen(request)
        txt = f.read()
        result = simplejson.loads(txt)
        if result.get('status') == 'error':
            print("Login error")
        session_string = result["session"]

        # Upload the text file to request a WCS solution
        upload_url = apiurl + "upload"
        upload_args = {
            'allow_commercial_use': 'd',
            'allow_modifications': 'd',
            'publicly_visible': 'y',
            'scale_units': 'arcsecperpix',
            'scale_type': 'ul',
            'scale_lower': self.pixel_scale + 0.05,       # arcsec/pix
            'scale_upper': self.pixel_scale - 0.05,       # arcsec/pix
            'parity': 0,
            'session': session_string
        }
        upload_json = simplejson.dumps(upload_args)

        f = open(upload_filename, 'rb')
        file_args = (upload_filename, f.read())

        m1 = MIMEBase('text', 'plain')
        m1.add_header('Content-disposition', 'form-data; name="request-json"')
        m1.set_payload(upload_json)

        m2 = MIMEApplication(file_args[1], 'octet-stream', encode_noop)
        m2.add_header('Content-disposition',
                      'form-data; name="file"; filename="%s"' % file_args[0])

        mp = MIMEMultipart('form-data', None, [m1, m2])

        fp = StringIO()
        g = MyGenerator(fp)
        g.flatten(mp)
        upload_data = fp.getvalue()
        upload_headers = {'Content-type': mp.get('Content-type')}

        # if False:
        #     print 'Sending headers:'
        #     print ' ', headers
        #     print 'Sending data:'
        #     print data[:1024].replace('\n', '\\n\n').replace('\r', '\\r')
        #     if len(data) > 1024:
        #         print '...'
        #         print data[-256:].replace('\n', '\\n\n').replace('\r', '\\r')
        #         print

        request = Request(url=upload_url, headers=upload_headers, data=upload_data)

        f = urlopen(request)
        txt = f.read()
        result = simplejson.loads(txt)
        stat = result.get('status')
        if stat == 'error':
            print("Upload error")
        submission_int = result["subid"]

        time.sleep(5)

        # Check submission status
        subcheck_url = apiurl + "submissions/" + str(submission_int)
        request = Request(url=subcheck_url)
        still_processing = True
        n_failed_attempts = 0
        while still_processing and n_failed_attempts < 5:
            try:
                f = urlopen(request)
                txt = f.read()
                result = simplejson.loads(txt)
                # print result
                still_processing = False
            except Exception:
                print("Submission doesn't exist yet, sleeping for 5s.")
                time.sleep(5)
                n_failed_attempts += 1
        if n_failed_attempts > 5:
            print("The submitted job has apparently timed out")

        job_id_list = result["jobs"]
        n_jobs = len(job_id_list)
        time.sleep(5)

        still_processing = True
        n_failed_attempts = 0
        n_failed_jobs = 0

        while still_processing and n_failed_attempts < 5 and n_failed_jobs < n_jobs:
            time.sleep(5)
            for job_id in job_id_list:
                jobcheck_url = apiurl + "jobs/" + str(job_id)
                request = Request(url=jobcheck_url)
                f = urlopen(request)
                txt = f.read()
                result = simplejson.loads(txt)
                print("Checking astrometry.net job ID"), job_id, result
                if result["status"] == "failure":
                    n_failed_jobs += 1
                    job_id_list.remove(job_id)
                if result["status"] == "success":
                    solved_job_id = job_id
                    still_processing = False
                    print("SOLVED")
            n_failed_attempts += 1

        if still_processing:
            print("Astrometry.net took too long to process.")

        if not still_processing:
            # Download WCS file
            os.system("wget -q  --output-document=" + output_path + "B1_wcs.fits" + " http://nova.astrometry.net/wcs_file/" + str(solved_job_id))
            os.system("wget -q  --output-document=" + output_path + "B1_corr.fits" + " http://nova.astrometry.net/corr_file/" + str(solved_job_id))
