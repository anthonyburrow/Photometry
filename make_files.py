import os.path

class MakeFiles:

	def __init__(self, root):
		self.root = root

	def ObsList():
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