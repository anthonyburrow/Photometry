import os.path


class MakeFiles:

    def __init__(self, root):
        self.root = root

    def StandsList(self):
        print("Listing photometry files for standard stars...")

        if not os.path.exists(self.root):
            os.makedirs(self.root)

        if os.listdir(self.root) != []:
            for date in sorted(os.listdir(self.root)):
                if os.path.isdir(os.path.join(self.root, date)) and date[:3] == "201":
                    G = open(self.root + date + "/mag_files.txt", 'w')

                    if os.listdir(self.root + date + "/") != []:
                        for file in sorted(os.listdir(self.root + date + "/")):
                            if file[-6:] == ".mag.1":
                                G.write(file + "\n")

                    G.close()
        else:
            print("No data in photometry directory.")
