import os.path

class MakeFiles:
    """Creates a MakeFiles object, controlling any files to be made.

    Creates any miscellaneous files to document observation and photometry information.

    Attributes:
        root: Root directory of photometry.
    """

	def __init__(self, root = "../photometry/"):
		self.root = root

	def ObsList():
	    """Creates observation lists.

        Creates a list of clusters observed and for each of these a list of observation
        dates which have completed photometry.

        """
		print ("Listing observations with complete photometry...")

		F = open(root + "obs_clusters.txt", 'w')

		for cluster in sorted(os.listdir(root)):
			if isdir(join(root, cluster)) and cluster[:3] == "NGC":
				F.write(cluster + "\n")

				G = open(root + cluster + "obs_clusters.txt", 'w') # May need an extra "/" before obs_clusters.txt

				for date in sorted(os.listdir(root + cluster + "/")):
					if isdir(join(root + cluster + "/", date)) and date[:3] == "201" and
					   isfile(root + cluster + "/" + date + "/B1.als.1") and
					   isfile(root + cluster + "/" + date + "/B3.als.1") and
					   isfile(root + cluster + "/" + date + "/V1.als.1") and
					   isfile(root + cluster + "/" + date + "/V3.als.1") and
					   isfile(root + cluster + "/" + date + "/R1.als.1") and
					   isfile(root + cluster + "/" + date + "/R3.als.1") and
					   isfile(root + cluster + "/" + date + "/H1.als.1") and
					   isfile(root + cluster + "/" + date + "/H3.als.1"):
						G.write(date + "\n")

				G.close()

		F.close()