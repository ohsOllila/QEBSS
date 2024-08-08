#!/usr/bin/python3

import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

if len(sys.argv) > 2:
    trjfile = sys.argv[1]
    topfile = sys.argv[2]
    trj = md.load(trjfile, top = topfile)
else:
    trjfile = sys.argv[1]
    trj = md.load(trjfile)

x, y = md.compute_contacts(trj)
z = md.geometry.squareform(x,y)

fig, ax = plt.subplots()
cmap = plt.cm.viridis
plot = ax.imshow(np.mean(z, axis = 0)[:,::-1], cmap = cmap)
ax.set_xlabel("Residue")
ax.set_ylabel("Residue")
cbar = plt.colorbar(plot, ax = ax)
cbar.set_label("Distance (nm)")

plt.show()
plt.savefig('mdmat.png')
