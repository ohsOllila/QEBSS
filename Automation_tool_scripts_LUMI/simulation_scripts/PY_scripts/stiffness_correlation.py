#!/usr/bin/python3

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob
import csv



def get_unit_vector(atom1, atom2):
    vector = atom2.position - atom1.position
    return vector / np.linalg.norm(vector)

def exponential_decay(s, A, k):
    return A * np.exp(-s / k)

def calculate_orientational_correlation(groFILE, trjFILE, frame_interval=1000):
    u = mda.Universe(groFILE, trjFILE)
    CAatoms = u.select_atoms("name CA")

    num_atoms = len(CAatoms.indices)
    max_s = num_atoms - 1

    orientational_correlations = np.zeros(max_s)
    counts = np.zeros(max_s)

    for frame in range(0, len(u.trajectory), frame_interval):
        u.trajectory[frame]
        for i in range(0, max_s):
            for s in range(0, max_s - i):
                u1 = get_unit_vector(CAatoms[i], CAatoms[i + 1])
                u2 = get_unit_vector(CAatoms[i + s], CAatoms[i + s + 1])

                orientational_correlations[s] += np.dot(u1, u2)
                counts[s] += 1

    avg_orientational_correlation = orientational_correlations / counts
    return avg_orientational_correlation

def plot_orientational_correlation(files):
    plt.figure(figsize=(10, 6))
    color_list = plt.cm.tab10.colors
    num_colors = len(color_list)
    persistence_data = []

   
    for groFILE, trjFILE in files:
        u = mda.Universe(groFILE, trjFILE)
        CAatoms = u.select_atoms("name CA")
        protein_length = len(CAatoms) - 1
        avg_correlation = calculate_orientational_correlation(groFILE, trjFILE)

        s_values_short = np.arange(0, 5)
        avg_correlation_short = avg_correlation[:5]
        params, _ = curve_fit(exponential_decay, s_values_short, avg_correlation_short, p0=[1, 1])

        persistence_data.append((params[1], groFILE, trjFILE, avg_correlation, params))

 
    persistence_data.sort(reverse=True, key=lambda x: x[0])

   
    for idx, (k, groFILE, trjFILE, avg_correlation, params) in enumerate(persistence_data):
        name_lab = '/'.join(trjFILE.split("/")[-4:-1])
        print(name_lab)
        for value in avg_correlation:
            print(value)

        color = color_list[idx % num_colors]  
        plt.plot(
            np.arange(0, len(avg_correlation)),
            avg_correlation, 
            lw=3,
            linestyle='-', 
            color=color, 
            label=f'{name_lab} Data (k={k:.2f})'
        )
        s_values_short = np.arange(0, 5)
        plt.plot(
           s_values_short,
           exponential_decay(s_values_short, *params),
           linestyle='--',
           lw=3,
           color=color
        )

#    plt.xscale('log')
#    plt.yscale('log')
    plt.xlabel("Residue Separation $s$")
 #   plt.xlim(-2, 25)
    plt.ylabel(r"Orientational Correlation $C(s)$")
    plt.title("Orientational Correlation Functions for Multiple Trajectories")
    plt.legend(fontsize=8, loc='best')
    plt.grid(True)
    plt.tight_layout()

    plt.savefig('orientational_correlation_sorted.png', dpi=600)
    plt.show()

path = "/scratch/project_2003809/cmcajsa/forcefield/"
selected_sim = [
    "aSyn/replica_01/AMBER03WS/",
    "ICL2/replica_05/DESAMBER/",
    "ChiZ/replica_04/AMBER99SB-DISP",
    "SNARE/replica_02/AMBER99SBWS",
    "SNARE/replica_03/AMBER99SBWS",
    "SNARE/replica_05/AMBER99SBWS",
    "KRS/replica_02/AMBER99SBWS"
]

if __name__ == "__main__":
    gro_files = []
    trj_files = []

    for sim in selected_sim:
        gro_files.extend(glob.glob(f"{path}{sim}/temp*gro"))  # Replace with the correct path and pattern
        trj_files.extend(glob.glob(f"{path}{sim}/md*noPBC*xtc"))
    files = list(zip(gro_files, trj_files))

    plot_orientational_correlation(files)
