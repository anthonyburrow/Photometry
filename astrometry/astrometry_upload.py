import sys
import time
import json
from urllib.parse import urlencode
from urllib.request import urlopen, Request, urlretrieve
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import os.path


class AstrometryUpload:

    def __init__(self):
        self.apiurl = 'http://nova.astrometry.net/api/'
        self.arcsecperpix = 0.465   # arcsec/pix
        self.cluster = 'NGC7419'

        self.AstrometryLogIn()

        self.fail_log = open('data/' + self.cluster + '_failures.txt', 'w')

    def AstrometryLogIn(self):
        # Log in to service
        print("Logging in...")
        with open('etc/astrometry_apikey') as G:
            apikey = G.readline()
        login_url = self.apiurl + 'login'
        data = {'request-json': json.dumps({"apikey": apikey})}
        data = urlencode(data).encode('utf-8')
        R = Request(login_url, data=data)
        f = urlopen(R)
        txt = f.read()
        result = json.loads(txt)
        stat = result.get('status')
        if stat == 'error':
            print("Error logging in.\n")
            sys.exit()
        else:
            print("Log in successful.\n")

        self.session_string = result["session"]

    def Process(self):
        filename = 'data/' + self.cluster + '_id.txt'
        data = np.genfromtxt(filename, dtype='str')

        for item in data:
            path = '../photometry/' + self.cluster + '/' + item[0] + '/'
            img = item[1]
            ID = item[2]
            wcs_filename = path + img + '_wcs.fits'
            corr_filename = path + img + '_corr.fits'
            if not os.path.isfile(wcs_filename) or not os.path.isfile(corr_filename):
                self.UploadFiles(path=path, img=img, ID=ID)

    def UploadFiles(self, path, img, ID):
        # Get default values for faster processing
        try:
            F = fits.getheader(path + img + '.fits')
            binning = F['XBINNING']
            centerRA = F['RA']
            centerDEC = F['DEC']
            # Reformat
            coo = SkyCoord(centerRA + centerDEC, unit=(u.hourangle, u.deg))
            centerRA = coo.ra.deg
            centerDEC = coo.dec.deg
        except Exception:
            binning = 1
            centerRA = None
            centerDEC = None

        # Submit url for processing
        upload_url = self.apiurl + "url_upload"
        upload_args = {
            "session": self.session_string,
            "url": 'https://drive.google.com/uc?export=download&id=' + ID,
            'publicly_visible': 'n',
            "scale_units": "arcsecperpix",
            "scale_lower": binning * self.arcsecperpix - 0.2,
            "scale_upper": binning * self.arcsecperpix + 0.2,
        }
        if centerRA is not None and centerDEC is not None:
            upload_args["center_ra"] = centerRA
            upload_args["center_dec"] = centerDEC
            upload_args["radius"] = 1

        data = {'request-json': json.dumps(upload_args)}
        data = urlencode(data).encode('utf-8')

        R = Request(upload_url, data=data)
        f = urlopen(R)
        txt = f.read()
        result = json.loads(txt)
        stat = result.get('status')
        if stat == 'error':
            print("Upload error for " + path + img)
            self.failed.append(path + img)
            return
        else:
            print("Upload started.")
        submission_int = result["subid"]

        time.sleep(5)

        # Check submission status
        subcheck_url = self.apiurl + "submissions/" + str(submission_int)
        request = Request(url=subcheck_url)
        still_processing = True
        n_failed_attempts = 0
        n_jobs = 0
        while n_failed_attempts < 30 and still_processing:
            try:
                f = urlopen(request)
                txt = f.read()
                result = json.loads(txt)
                n_jobs = len(result["jobs"])
                # print("Result: \n", result)
                if n_jobs > 0 and result["jobs"] != [None]:
                    still_processing = False
                    print("Upload finished...")
            except Exception:
                print("Submission doesn't exist yet, sleeping for 5s.")
                n_failed_attempts += 1
            time.sleep(5)

        if n_failed_attempts > 5:
            print("The submitted job has timed out for " + path + img)
            self.failed.append(path + img)
            return

        print("Processing...")

        job_id_list = result["jobs"]
        # print(result)
        n_jobs = len(job_id_list)

        still_processing = True
        n_failed_attempts = 0
        n_failed_jobs = 0

        while still_processing and n_failed_attempts < 30:   # and n_failed_jobs < n_jobs:
            time.sleep(5)
            for job_id in job_id_list:
                jobcheck_url = self.apiurl + "jobs/" + str(job_id)
                request = Request(url=jobcheck_url)
                f = urlopen(request)
                txt = f.read()
                # print(txt)
                result = json.loads(txt)
                # print(result)
                # print("Checking astrometry.net job ID", job_id, result)
                if result["status"] == "failure":
                    n_failed_jobs += 1
                    job_id_list.remove(job_id)
                if result["status"] == "success":
                    solved_job_id = job_id
                    still_processing = False
                    # print("SOLVED")
            n_failed_attempts += 1

        if still_processing:
            print("Astrometry.net took too long to process for " + path + img + '\n')
            self.fail_log.write(path + img + ' ' + ID + '\n')
            return
        else:
            urlretrieve("http://nova.astrometry.net/wcs_file/" + str(solved_job_id), path + img + '_wcs.fits')
            urlretrieve("http://nova.astrometry.net/corr_file/" + str(solved_job_id), path + img + '_corr.fits')
            print(img + " files saved in " + path + "\n")


if __name__ == "__main__":
    process = AstrometryUpload()
    process.Process()
