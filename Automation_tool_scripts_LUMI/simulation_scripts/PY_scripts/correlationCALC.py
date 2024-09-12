#!/usr/bin/python3


import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import glob

groFILE = glob.glob('temp_*.gro')[0]
trjFILE = glob.glob('md*noPBC*xtc')[0]
outfile_png = 'LRAEcorrelationHELICALstart.png'
outfile_csv = 'LRAEcorrelationHELICALstart.csv'

u = mda.Universe(groFILE, trjFILE)
CAatoms = u.select_atoms("name CA")


num_atoms = len(CAatoms.ix)
num_vectors = num_atoms - 1
vec = np.zeros((num_vectors, 3))
matrix = np.zeros((num_vectors, num_vectors))

for frame in u.trajectory:
	for i in range(num_vectors):
		vec[i]=CAatoms[i+1].position-CAatoms[i].position
		vec[i] /= np.linalg.norm(vec[i])

	matrix += np.dot(vec, vec.T)

matrix /= len(u.trajectory)


np.savetxt(outfile_csv, matrix, delimiter=',')


color_map = plt.imshow(matrix,vmin=-1, vmax=1, cmap="seismic", origin='lower')

cbar = plt.colorbar(color_map)
cbar.ax.tick_params(labelsize=16)

plt.tick_params(axis='both', which='major', labelsize=12)
plt.xlabel(r"C$_{\alpha}$ Pair Index", fontsize=16)
plt.ylabel(r"C$_{\alpha}$ Pair Index", fontsize=16)



plt.tight_layout()

plt.savefig(outfile_png, dpi=600)
plt.close("all")




