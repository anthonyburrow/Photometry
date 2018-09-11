import numpy as np
from zeropoint import ZeroPoint
import os.path


class ApertureCorrect:

    def __init__(self, cluster, outputDate, date, totalAperture, varAperture):
        self.cluster = cluster
        self.outputDate = outputDate
        self.date = date
        self.totalAperture = totalAperture
        self.varAperture = varAperture

    def GetCorrection(self):
        # Read standard star data for each aperture
        files = ZeroPoint(self.date, self.totalAperture).StandsList()
        Bfiles = [file for file in files if file[0] == "B"]
        Vfiles = [file for file in files if file[0] == "V"]
        Rfiles = [file for file in files if file[0] == "R"]
        files = [Bfiles, Vfiles, Rfiles]

        BSample = []
        VSample = []
        RSample = []

        for xFiles in files:
            for file in xFiles:
                filename = "../../standards/" + self.date + "/ap" + str(self.totalAperture) + "/" + file
                totalData = ZeroPoint(self.date, self.totalAperture).magRead(filename)
                filename = "../../standards/" + self.date + "/ap" + str(self.varAperture) + "/" + file
                varData = ZeroPoint(self.date, self.varAperture).magRead(filename)

                for i in range(0, len(totalData)):
                    # To correct from bigger aperture to smaller, offset is (varAperture - totalAperture)
                    diff = float(varData[i][2]) - float(totalData[i][2])

                    if file[0] == 'B':
                        BSample.append(diff)
                    elif file[0] == 'V':
                        VSample.append(diff)
                    elif file[0] == 'R':
                        RSample.append(diff)

        BApertureCorrection = np.mean(BSample)
        VApertureCorrection = np.mean(VSample)
        RApertureCorrection = np.mean(RSample)

        data = [BApertureCorrection, VApertureCorrection, RApertureCorrection]

        # Write to file
        if not os.path.exists("../../standards/" + self.outputDate + "/"):
            os.makedirs("../../standards/" + self.outputDate + "/")

        filename = "../../standards/" + self.outputDate + "/" + self.cluster + "_aperture_corrections.dat"
        with open(filename, 'w') as F:
            np.savetxt(F, data, fmt='%.3f')


if __name__ == "__main__":
    cluster = 'NGC663'
    obs = [
        ["20150829", "20151102", 20.0, 6.0],
        ["20151102", "20151102", 20.0, 6.0],
        ["20151104", "20151104", 20.0, 6.5],
        ["20151106", "20151106", 16.0, 5.0],
        ["20151201", "20151106", 16.0, 5.0],
        ["20151204", "20151106", 16.0, 5.0],
        ["20161017", "20161018", 20.0, 5.5],
        ["20161018", "20161018", 20.0, 5.5],
        ["20161019", "20161019", 16.0, 6.5],
        ["20161021", "20161021", 20.0, 8.5],
        ["20161023", "20161021", 20.0, 8.5],
        ["20161213", "20161021", 20.0, 8.5]
    ]

    for date in obs:
        ApertureCorrect(cluster, date[0], date[1], date[2], date[3]).GetCorrection()

    cluster = 'NGC869'
    obs = [
        ["20150829", "20151102", 20.0, 6.5],
        ["20151102", "20151102", 20.0, 6.5],
        ["20151104", "20151104", 20.0, 8.0],
        ["20151105", "20151104", 20.0, 8.0],
        ["20151106", "20151106", 16.0, 6.0],
        ["20151201", "20151106", 16.0, 6.0],
        ["20151204", "20151106", 16.0, 6.0],
        ["20161017", "20161018", 20.0, 5.0],
        ["20161018", "20161018", 20.0, 5.0],
        ["20161019", "20161019", 16.0, 8.0],
        ["20161021", "20161021", 20.0, 6.0],
        ["20161023", "20161021", 20.0, 6.0]
    ]

    for date in obs:
        ApertureCorrect(cluster, date[0], date[1], date[2], date[3]).GetCorrection()

    cluster = 'NGC884'
    obs = [
        ["20150827", "20151102", 20.0, 6.5],
        ["20150829", "20151102", 20.0, 6.5],
        ["20150830", "20151102", 20.0, 6.5],
        ["20151102", "20151102", 20.0, 6.5],
        ["20151104", "20151104", 20.0, 6.0],
        ["20151106", "20151106", 16.0, 4.5],
        ["20151201", "20151106", 16.0, 4.5],
        ["20151204", "20151106", 16.0, 4.5],
        ["20161017", "20161018", 20.0, 5.0],
        ["20161018", "20161018", 20.0, 5.0],
        ["20161019", "20161019", 16.0, 12.0],
        ["20161021", "20161021", 20.0, 8.0],
        ["20161023", "20161021", 20.0, 8.0]
    ]

    for date in obs:
        ApertureCorrect(cluster, date[0], date[1], date[2], date[3]).GetCorrection()
