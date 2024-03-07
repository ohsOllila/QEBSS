#!/usr/bin/python3

import glob
import os
from pymol import cmd
import matplotlib.pyplot as plt
import numpy as np



md=glob.glob("*noPBC*xtc*")
temp=glob.glob("temp*gro")

cmd.load(temp[0])
cmd.load_traj(md[0])
cmd.do("remove ss h")
cmd.do("remove first name CA")
cmd.do("select first_atom,  first name CA")
cmd.do("remove last name CA")


list=["008", "*010", "012", "014", "018", "031", "037", "038", "044", "047", "050","051", "053", "066", "072", "075", "096", "100", "101", "104", "105", "106","109", "110", "111", "114", "118", "122", "126", "128", "131", "147", "150", "155", "158", "159", "160", "164", "169", "171", "173"]


for i in list:
	if i in md[0]:
		cmd.do("remove last name CA")


if "080" in md[0]:
	cmd.do("remove last name CA")
	cmd.do("remove last name CA")


		
cmd.do("select last_atom, last name CA")

distances=[]
for a in range(1, 10002):
	st0=cmd.distance("(last_atom)", "(first_atom)", state=a)
	distances.append(st0)
	
directory_path = os.getcwd()
folder_name = os.path.basename(directory_path)

f=open(folder_name + '_distance_variation_in_helix_connected_residue.txt','w')	
for items in distances:
	f.write("%8.3f\n"%items)
f.close()

x=[]
for number in np.arange(0, 100.1, 0.1 ):
    x.append(number)
    
y=distances[0:10002:10]


plt.plot(x, y)

plt.xlabel('Time, ns')
plt.ylabel('Distance, Å')
plt.margins(x=0)
plt.title("Distance between residues connected to helices over time")
plt.savefig(folder_name + '_distance_variation_in_helix_connected_residue.png')
