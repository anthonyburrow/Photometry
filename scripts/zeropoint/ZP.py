import numpy

#read in list of files

files = numpy.loadtxt("tempfiles/files.txt",dtype = 'string')

#read in data

def navigate(x):
    return 'tempfiles/' + x

class Data:
    
    def __init__(self, dataset):

        self.xcoo  = dataset[:,0]
        self.ycoo  = dataset[:,1]
        self.mag   = dataset[:,2]
        self.merr  = dataset[:,3]

fileno = len(files)
bfiles = [elem for elem in files if elem[0] == "b"]
rfiles = [elem for elem in files if elem[0] == "r"]
vfiles = [elem for elem in files if elem[0] == "v"]

for i in range(0,fileno):

    data = files[i]
    z = numpy.genfromtxt(navigate(data), usecols = (1,2,33,34), dtype = 'f')
    x = Data(z)
    exec(data + " = x")

#match between files

tol = 10
counter = 0
starno = len(b1.xcoo)

#print b1.mag
#print b2.mag
#print r1.mag
#print r2.mag
#print v1.mag
#print v2.mag

Bmean = []
Rmean = []
Vmean = []

for i in range (0,starno):

    x0 = str(b1.xcoo[i])
    y0 = str(b1.ycoo[i])  

    for j in range (0,starno):

        for k in range (1,fileno):

            c = files[k]

            x1 = c + ".xcoo[" + str(j) + "]" 
            y1 = c + ".ycoo[" + str(j) + "]" 

            exec("xmag = abs(" + x0 + "-" + x1 + ")")
            exec("ymag = abs(" + y0 + "-" + y1 + ")")

            if xmag <= tol and ymag <= tol:
                                            
                if c[0] == "r":
        
                    Rmag = []
                    
                    for elem in rfiles:

                        exec("z = "+ elem + ".mag[" + str(j) + "]") 

                        if len(Rmag) > 0:

                            m = numpy.mean(Rmag)

                        if len(Rmag) == 0 or abs(z - m) <= 0.4:
                        
                            exec("Rmag.append(" + str(z)  + ")")

                        else:

                            Rmag = []

                    if len(Rmag) != 0:        

                        m = numpy.mean(Rmag)

                    if m not in Rmean and len(Rmag) != 0:
                        
                        Rmean.append(m)

                    #print Rmag

                elif c[0] == "b":
        
                    Bmag = []
                    
                    for elem in bfiles:

                        exec("z = "+ elem + ".mag[" + str(j) + "]") 

                        if len(Bmag) > 0:

                            m = numpy.mean(Bmag)

                        if len(Bmag) == 0 or abs(z - m) <= 0.4:
                        
                            exec("Bmag.append(" + str(z)  + ")")

                        else:

                            Bmag = []

                    if len(Bmag) != 0:        

                        m = numpy.mean(Bmag)

                    if m not in Bmean and len(Bmag) != 0:
                        
                        Bmean.append(m)

                    print Bmag

#print Rmean
print Bmean
print Rmean
