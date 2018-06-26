import numpy as np


class LowError:

    def __init__(self, cluster, date, app, maxError=0.04):
        self.cluster = cluster
        self.date = date
        self.app = app
        self.maxError = maxError

    def Process(self, data):
        lowError_data = []

        Rerr = np.array(data)[:, 7]
        Herr = np.array(data)[:, 9]

        R_Herr = np.sqrt(Rerr**2 + Herr**2)

        for i in range(0, len(data)):
            if R_Herr[i] < self.maxError:
                lowError_data.append(data[i])

        # Output to file
        filename = "../output/" + self.cluster + "/" + self.date + "/phot_" + self.app.phot_type + "_lowError.dat"
        with open(filename, 'w') as F:
            for item in data:
                for value in item:
                    F.write(str(value) + " ")
                F.write("\n")

        return lowError_data
