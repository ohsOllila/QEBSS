### plot distancemaps

import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

base_dir = os.path.dirname(os.getcwd()) + "/"
RESULTS = base_dir + 'results/'


proteins=["SNARE", "ChiZ", "ICL2", "aSyn", "KRS"]
models = ["model_01", "model_02", "model_03", "model_04", "model_05"]
forcefields = ["AMBER03WS", "AMBER99SB-DISP", "AMBER99SBWS", "CHARMM36M", "DESAMBER"]


fig, axs = plt.subplots(25, 5, figsize=(40, 40))
for k, prot in enumerate(proteins):
    for i, model in enumerate(models):
        for j, forcefield in enumerate(forcefields):
            try:
                file_path = base_dir + prot + "/" + model + "/" + forcefield + "/LRAEcorrelationHELICALstart.csv"
                df = pd.read_csv(file_path, index_col=0)
                axs[5*k+i, j].imshow(df.values, vmin=-1, vmax=1, cmap="seismic", origin='lower')
            except:
                axs[5*k + i, j].imshow(np.zeros((10, 10)), cmap='gray')

#plt.tight_layout()
fig.subplots_adjust(right=0.9, top=0.9)

fig.savefig(RESULTS + "major_correlation.png", dpi=300, bbox_inches='tight')
plt.close(fig)


