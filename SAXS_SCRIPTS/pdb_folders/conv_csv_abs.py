import pandas as pd
import numpy as np

# Step 1: Load the CSV file
input_csv = "exp.csv"  # Replace with your input CSV file name
output_abs = "exp_75mM.abs"  # Replace with your desired ABS file name

# Load the CSV
df = pd.read_csv(input_csv)

df.iloc[:, 0] = df.iloc[:, 0] / 10


# Step 3: Save to ABS format
# Define the ABS format (e.g., space-separated with no headers or index)
df.to_csv(output_abs, index=False, header=False, sep=' ')
