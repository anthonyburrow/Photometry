#This program reorganizes .mag files into one-line format, from a directory to itself.
#Currently you must rename all input files and put as many iterations as there are files,
#  and rename all variables.


B1 = open('B_015355.mag.1')
b1_band = B1.readlines()[75:]
B1.close()

B2 = open('B_020056.mag.1')
b2_band = B2.readlines()[75:]
B2.close()

R1 = open('R_022335.mag.1')
r1_band = R1.readlines()[75:]
R1.close()

R2 = open('R_022940.mag.1')
r2_band = R2.readlines()[75:]
R2.close()

V1 = open('V_020745.mag.1')
v1_band = V1.readlines()[75:]
V1.close()

V2 = open('V_021342.mag.1')
v2_band = V2.readlines()[75:]
V2.close()

b1_data = []
b2_data = []
r1_data = []
r2_data = []
v1_data = []
v2_data = []

for i in range(0,int(len(b1_band)/5.)):
	list = b1_band[5*i] + ' ' + b1_band[5*i+1] + ' ' + b1_band[5*i+2] + ' ' + b1_band[5*i+3] + ' ' + b1_band[5*i+4]
	b1_data.append(list)

for i in range(0,int(len(b2_band)/5.)):
	list = b2_band[5*i] + ' ' + b2_band[5*i+1] + ' ' + b2_band[5*i+2] + ' ' + b2_band[5*i+3] + ' ' + b2_band[5*i+4]
	b2_data.append(list)

for i in range(0,int(len(r1_band)/5.)):
	list = r1_band[5*i] + ' ' + r1_band[5*i+1] + ' ' + r1_band[5*i+2] + ' ' + r1_band[5*i+3] + ' ' + r1_band[5*i+4]
	r1_data.append(list)

for i in range(0,int(len(r2_band)/5.)):
	list = r2_band[5*i] + ' ' + r2_band[5*i+1] + ' ' + r2_band[5*i+2] + ' ' + r2_band[5*i+3] + ' ' + r2_band[5*i+4]
	r2_data.append(list)

for i in range(0,int(len(v1_band)/5.)):
	list = v1_band[5*i] + ' ' + v1_band[5*i+1] + ' ' + v1_band[5*i+2] + ' ' + v1_band[5*i+3] + ' ' + v1_band[5*i+4]
	v1_data.append(list)

for i in range(0,int(len(v2_band)/5.)):
	list = v2_band[5*i] + ' ' + v2_band[5*i+1] + ' ' + v2_band[5*i+2] + ' ' + v2_band[5*i+3] + ' ' + v2_band[5*i+4]
	v2_data.append(list)


F = open('tempfiles/b1','w')

for data_b in b1_data:
	curr_b1 = data_b.split()
        F.write(' '.join(curr_b1))
        F.write('\n')

F.close()

F = open('tempfiles/b2','w')

for data_b in b2_data:
       curr_b2 = data_b.split()
       F.write(' '.join(curr_b2))
       F.write('\n')

F.close()

F = open('tempfiles/r1','w')

for data_r in r1_data:
	curr_r1 = data_r.split()
        F.write(' '.join(curr_r1))
        F.write('\n')

F.close()

F = open('tempfiles/r2','w')

for data_r in r2_data:
	curr_r2 = data_r.split()
        F.write(' '.join(curr_r2))
        F.write('\n')

F.close()

F = open('tempfiles/v1','w')

for data_v in v1_data:
	curr_v1 = data_v.split()
        F.write(' '.join(curr_v1))
        F.write('\n')

F.close()

F = open('tempfiles/v2','w')

for data_v in v2_data:
	curr_v2 = data_v.split()
        F.write(' '.join(curr_v2))
        F.write('\n')

F.close()
