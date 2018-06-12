import MakeFiles
import Match
import BeFilter
import Plot

root = "../photometry/"

makeFiles = MakeFiles(root)
makeFiles.ObsList()

autoCheck = raw_input("Process every night? [y/n]: ") or "y"

if autoCheck = "y":
	with open("../photometry/obs_clusters.txt") as F:
	    for cluster in F:
	        with open("../photometry/" + cluster + "/obs_dates.txt") as G:
        		for date in G:
        			ProcessDate(cluster, str(date))
        			
        	ProcessCluster(cluster)

elif autoCheck = "n":
	cluster = raw_input("Cluster: ")
	date = raw_input("Date: ")
	ProcessDate(cluster, date)

def ProcessDate(cluster, date):
	input_directory = "../photometry/" + cluster + "/" + str(date) + "/"
	output_directory = "../output/" + cluster + "/" + str(date) + "/"

    print ("Compiling all data for " cluster + " on " + date + "...")

	match = Match(output_directory, input_directory, "psf", 5.0, 0.5)
	match.LowError()

	print ("Extracting Be candidate data...")

	beFilter = BeFilter(output_directory, "psf", true)
	beFilter.LowError()

	print ("Generating plots...")

	plot = Plot(output_directory, true, true)
	plot.ColorMagnitudeDiagram()
	plot.TwoColorDiagram()

def ProcessCluster(cluster):
	print("Scaling observation nights for " + cluster)

	scale = Scale(cluster, "psf", 10)
