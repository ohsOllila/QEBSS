### plot distancemaps

import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

base_dir = os.path.dirname(os.getcwd()) + "/"


models = ["model_01", "model_02", "model_03", "model_04", "model_05"]
forcefields = ["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]


# location of CSV with average distances:
# base_dir/system/model/forcefield/CA_avg_distances.csv

# Function to create the plot
def create_avg_distances_plot(system):
    fig, axs = plt.subplots(5, 5, figsize=(20, 20))

    #fig.suptitle(f"{system} - Average CA Distances", fontsize=24, y=0.95)

    # Create a list to store all imshow objects and all data
    ims = []
    all_data = []

    # First pass: collect all data to determine global vmin and vmax
    for i, model in enumerate(models):
        for j, forcefield in enumerate(forcefields):
            # Construct the path to the CSV file
            file_path = glob.glob(base_dir + "ChiZ/" + model + "/" + forcefield + "/CA_avg_distances.csv")[0]
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path, index_col=0)
            all_data.append(df.values)

    # Calculate global vmin and vmax
    vmin = np.floor(np.min([np.min(data) for data in all_data]) / 10) * 10
#    vmax = np.ceil(np.max([np.max(data) for data in all_data]) / 10) * 10
    vmax = 20
    # Second pass: create the heatmaps
    for i, model in enumerate(models):
        for j, forcefield in enumerate(forcefields):
            # Construct the path to the CSV file
            file_path = glob.glob(base_dir + "ChiZ/" + model + "/" + forcefield + "/CA_avg_distances.csv")[0]
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path, index_col=0)

            # Create the heatmap
            im = axs[i, j].imshow(df.values, cmap='jet', vmin=vmin, vmax=vmax)
            ims.append(im)

            # Set the title for each subplot
            axs[i, j].set_title(f"{model}\n{forcefield}", fontsize=18)

            # Set tick font size for all subplots
            axs[i, j].tick_params(axis='both', which='major', labelsize=14)

            # Add x-axis ticks to the bottom row
            if i == 4:
                axs[i, j].set_xticks(np.arange(0, len(df.columns), len(df.columns)//5))
                axs[i, j].set_xticklabels(np.arange(1, len(df.columns)+1, len(df.columns)//5))
            else:
                axs[i, j].set_xticks([])

            # Add y-axis ticks to the leftmost column
            if j == 0:
                axs[i, j].set_yticks(np.arange(0, len(df.index), len(df.index)//5))
                axs[i, j].set_yticklabels(np.arange(1, len(df.index)+1, len(df.index)//5))
            else:
                axs[i, j].set_yticks([])

            axs[i, j].invert_yaxis()

    # Adjust the layout and add space for the colorbar
    plt.tight_layout()
    fig.subplots_adjust(right=0.9, top=0.9)

    # Add a colorbar to the figure
    cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.3])
    cbar = fig.colorbar(ims[0], cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Distance (Å)', rotation=270, labelpad=20, fontsize=18)

    # Set colorbar ticks
    cbar.set_ticks(np.linspace(vmin, vmax, num=6))

    return fig, vmin, vmax


# Generate and save contact maps for each system
systems = glob.glob(base_dir + "ChiZ/")

# location of output PNG file for each system:
# base_dir/system/PLOTS/

for system in systems:
    fig, vmin, vmax = create_avg_distances_plot(system)
    path_parts=system.split("/")
    path_parts.insert(-2, "results")
    output_dir='/' + '/'.join(part.strip('/') for part in path_parts if part) + "/rep_to_exp_data"

    fig.savefig(os.path.join(output_dir, "Distance_map_combined.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    print("Figure saved in ", os.path.join(output_dir, "Distance_map_combined.png"))
    print(f"vmin: {vmin}, vmax: {vmax}")

