import os.path
import requests
import pandas as pd


class Zeropoint:

    def __init__(self, root, date):
        self.root = root
        self.date = date

    def StandsList(self):
        print("Listing photometry files for standard stars...")

        if not os.path.exists(self.root):
            os.makedirs(self.root)

        files = []
        if os.listdir(self.root) != []:
            if os.path.isdir(os.path.join(self.root, self.date)) and self.date[:3] == "201":
                if os.listdir(self.root + self.date + "/") != []:
                    for file in sorted(os.listdir(self.root + self.date + "/")):
                        if file[-6:] == ".mag.1":
                            files.append(file)
        else:
            print("There is no standard star photometry available.")

        if files == []:
            print("No standard star photometry for " + self.date)

        return files

    def magRead(self, filename):
        with open(self.root + self.date + "/" + filename) as F:
            file = F.readlines()[75:]

        data = []

        for i in range(0, int(len(file) / 5.)):
            # Concatenate lines
            combined = file[5 * i] + ' ' + file[5 * i + 1] + ' ' + file[5 * i + 2] + ' ' + file[5 * i + 3] + ' ' + file[5 * i + 4]
            # Select values needed in data set: X, Y, mag, mag error
            combined = combined.split()
            selected = []
            selected.extend((combined[7], combined[8], combined[33], combined[34]))
            data.append(selected)

        return data

    def GetStandardData():
        url = "http://james.as.arizona.edu/~psmith/61inch/ATLAS/charts/c121.html"
        html = requests.get(url).content
        df = pd.read_html(html)
        # ...

    def ZeroPoint(self, date):
        files = self.StandsList()
        Bfiles = [file for file in files if file[0] == "B"]
        Vfiles = [file for file in files if file[0] == "V"]
        Rfiles = [file for file in files if file[0] == "R"]
