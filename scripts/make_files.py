import os.path


class MakeFiles:
    """Creates a MakeFiles object, controlling any files to be made.

    Creates any miscellaneous files to document observation and photometry information.

    Attributes:
            root: Root directory of photometry.
    """

    def __init__(self, root):
        self.root = root

    def ObsList(self):
        """Creates observation lists.

        Creates a list of clusters observed and for each of these a list of observation
        dates which have completed photometry.

        """
        print("Listing observations with complete photometry...")

        if not os.path.exists(self.root):
            os.makedirs(self.root)

        F = open(self.root + "obs_clusters.txt", 'w')

        if os.listdir(self.root) != []:
            for cluster in sorted(os.listdir(self.root)):
                if os.path.isdir(os.path.join(self.root, cluster)) and cluster[:3] == "NGC":
                    F.write(cluster + "\n")

                    G = open(self.root + cluster + "/obs_clusters.txt", 'w')

                    if os.listdir(self.root + cluster + "/") != []:
                        for date in sorted(os.listdir(self.root + cluster + "/")):
                            if os.path.isdir(os.path.join(self.root + cluster + "/", date)) and date[:3] == "201" and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/B1.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/B3.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/V1.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/V3.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/R1.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/R3.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/H1.als.1") and \
                                    os.path.isfile(self.root + cluster + "/" + date + "/H3.als.1"):
                                G.write(date + "\n")

                    G.close()
        else:
            print("No data in photometry directory.")

        F.close()
