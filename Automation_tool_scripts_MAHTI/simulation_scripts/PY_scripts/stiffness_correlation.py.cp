#!/usr/bin/python3


import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob

def get_unit_vector(atom1, atom2):
    vector = atom2.position - atom1.position
    return vector / np.linalg.norm(vector)

def exponential_decay(s, A, lp):
    return A * np.exp(-s / lp)

def calculate_orientational_correlation(groFILE, trjFILE, frame_interval=1000):
    u = mda.Universe(groFILE, trjFILE)
    CAatoms = u.select_atoms("name CA")

    num_atoms = len(CAatoms.indices)
    max_s = num_atoms - 1

    orientational_correlations = np.zeros(max_s)
    counts = np.zeros(max_s)

    for frame in range(0, len(u.trajectory), frame_interval):
        u.trajectory[frame]  # This moves to the specified frame number
        for i in range(0, max_s):
            for s in range(1, max_s - i):
                u1 = get_unit_vector(CAatoms[i], CAatoms[i + 1])
                u2 = get_unit_vector(CAatoms[i + s], CAatoms[i + s + 1])

                orientational_correlations[s] += np.dot(u1, u2)
                counts[s] += 1

    avg_orientational_correlation = orientational_correlations / counts
    return avg_orientational_correlation

def plot_orientational_correlation(files):
    plt.figure(figsize=(10, 6))
    color_list=['red', 'blue', 'green', 'purple', 'orange']

    persistence_lengths = []
    labels = []

    for groFILE, trjFILE in files:
        u = mda.Universe(groFILE, trjFILE)
        CAatoms = u.select_atoms("name CA")  # Select C-alpha atoms
        protein_length = len(CAatoms) - 1
        avg_correlation = calculate_orientational_correlation(groFILE, trjFILE)

        s_values_short = np.arange(1, 25)
        avg_correlation_short = avg_correlation[1:25]
        params, _ = curve_fit(exponential_decay, s_values_short, avg_correlation_short, p0=[1, 1])

        persistence_lengths.append(params[1])
        labels.append(f'{trjFILE} (lp={params[1]:.2f} nm)')

        plt.plot(np.arange(1, len(avg_correlation[1:])+1), avg_correlation[1:], lw=4, linestyle='-', label=f'{trjFILE} Data')
#        plt.plot(s_values_short, exponential_decay(s_values_short, *params), '-', label=f'{trjFILE} (lp={params[1]:.2f} nm)')

#    sorted_indices = np.argsort(persistence_lengths)[::-1]
#    sorted_labels = [labels[i] for i in sorted_indices]
    
#    plt.xscale('log')  # Set x-axis to logarithmic scale
#    plt.yscale('log')
#    plt.xlim(1, 10)  # Cutoff at 10 for x-axis

    plt.xlabel("Residue Separation $s$")
    plt.ylabel(r"Orientational Correlation $C(s)$")
    plt.title("Orientational Correlation Functions for Multiple Trajectories")
    plt.legend(fontsize=5)
    plt.grid(True)

    plt.tight_layout()

    plt.savefig('orientational_correlation_multiple.png', dpi=600)
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
