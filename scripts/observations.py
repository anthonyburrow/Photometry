import os.path


class Observations:

    def __init__(self):
        self.root = "../photometry/"

        self.Check()

    def Check(self):
        if not os.path.exists(self.root):
            os.makedirs(self.root)

    def ListClusters(self):
        clusters = []
        if os.listdir(self.root) != []:
            for cluster in sorted(os.listdir(self.root)):
                if os.path.isdir(os.path.join(self.root, cluster)) and cluster[:3] == "NGC":
                    clusters.append(cluster)
        else:
            print("There are no files in the photometry directory.")

        if clusters == []:
            print("There are no appropriate cluster directories in the photometry directory.")

        return clusters

    def ListDates(self, cluster):
        dates = []
        if os.path.isdir(self.root + cluster + "/"):
            if os.listdir(self.root + cluster + "/") != []:
                for date in sorted(os.listdir(self.root + cluster + "/")):
                    files = ["B1.als.1", "B3.als.1", "V1.als.1", "V3.als.1", "R1.als.1", "R3.als.1", "H1.als.1", "H3.als.1"]
                    if os.path.isdir(os.path.join(self.root + cluster + "/", date)) and date[:3] == "201" and \
                            set(files).issubset(os.listdir(self.root + cluster + "/" + date + "/")):
                        dates.append(date)
            else:
                print(cluster + " directory does not have any files in it.")
        else:
            print(cluster + " is not a valid cluster directory.")

        return dates
