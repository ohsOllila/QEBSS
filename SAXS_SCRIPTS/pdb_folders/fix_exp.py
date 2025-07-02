import glob
import numpy as np
import os

# Folder containing your exp*.abs files
exp_files = glob.glob("exp*.abs")

for exp_file in exp_files:
    # Load the data
    data = np.loadtxt(exp_file, skiprows=1)
    
    # Multiply q-values (first column) by 10
    data[:, 0] *= 10

    # Read header to preserve it
    with open(exp_file, 'r') as f:
        header = f.readline()

    # Prepare output file path
    base_name = os.path.basename(exp_file)
    output_file = os.path.join(base_name.replace('.abs', '_fixed.abs'))

    # Save corrected data to new file
    np.savetxt(output_file, data, header=header.strip(), comments='', fmt="%.8e")
    print(f"Created corrected file: {output_file}")
