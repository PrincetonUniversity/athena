# python script to plot mesh structure in 'mesh_structure.dat' file produced
# by '-m #np' on command line

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x=[]
y=[]
z=[]
with open("mesh_structure.dat") as f:
    for line in f:
        if line[0]!='\n' and line[0]!='#':
                numbers_str = line.split()
                x.append(float(numbers_str[0]))
                y.append(float(numbers_str[1]))
		#[JMSHI
                #z.append(float(numbers_str[2]))
		if(len(numbers_str)>2):
                  z.append(float(numbers_str[2]))
		else:
		  z.append(0.0)
		#JMSHI]
        if line[0]=='\n' and len(x)!=0:
	        plt.plot(x,y,z,'k-')
		x=[]
		y=[]
		z=[]
plt.show()
