#!/usr/bin/env python

import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import glob
import shutil
from matplotlib.ticker import FixedLocator, FixedFormatter

def unquote(s):
	return s[1+s.find('"'):s.rfind('"')]

def uncomment(s):
	return s[2+s.find('/*'):s.rfind('*/')]
    
    
def col(c):
	color = c.split('/*')
	value = unquote(color[1])
	color = unquote(color[0]).split()
	return color[0], value

temp_file=glob.glob("*xpm")



xpm = open(temp_file[0])
png_file_name = temp_file[0].replace(".xpm", ".png")
csv_file_name = temp_file[0].replace(".xpm", ".csv")

	# Read in lines until we fidn the start of the array

meta = [xpm.readline()]

while not meta[-1].startswith("static char *gromacs_xpm[]"):
	meta.append(xpm.readline())
	    
	    
# The next line will contain the dimensions of the array

dim = xpm.readline()

# There are four integers surrounded by quotes

nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]

# The next dim[2] lines contain the color definitions

# Each pixel is encoded by dim[3] bytes, and a comment

# at the end of the line contains the corresponding value

colors = dict([col(xpm.readline()) for i in range(nc)])

mdmat = []

for i in xpm:
	if i.startswith("/*"):
		continue
	j = unquote(i)
	z = [float(colors[j[k:k+nb]]) for k in range(0,nx,nb)]
	mdmat.append(z)
	    
mdmat = np.array(mdmat)
np.savetxt(csv_file_name, mdmat[::-1,:], delimiter = ",")

fig, ax = plt.subplots()
cmap = plt.cm.viridis
plot = ax.imshow(mdmat[::-1,:], cmap = cmap, origin = "lower", vmin=0, vmax=1.4)
num_residues = np.shape(mdmat)[0]

ax.set_xlabel("Residue")
ax.set_ylabel("Residue")

cbar = plt.colorbar(plot, ax = ax)
cbar.set_label("Distance (nm)")
plt.savefig(png_file_name)
plt.close()

