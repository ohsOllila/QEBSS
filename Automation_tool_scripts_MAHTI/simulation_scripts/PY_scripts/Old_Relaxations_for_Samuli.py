#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""Jupyter Notebook for relaxation time analysis"""
#The main analysis parts adapted from script by  H. Antila, with help from S. Ollila and T. Ferreira
#saved in relaxation_times.py
# Last modified by R. Nencini, 19.10.2021
import yaml
import sys
import glob
import csv
import subprocess
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

#gyromagnetic ratios for further developmant
# !!! NOTICE!!!
#values taken from matlab code and projecct work and slightly different than those in Wikipedia
#these values are also in the external file --> if change is needed; has to be changed there
#values here in jupyter just for the information and verify, why they are different
#!!! NOTICE END !!!
gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;


# In[2]:


"""Parameters to be specified by the user"""
OP=0 # order parameter
smallest_corr_time=0 # enter in log scale -3 fs; 0 ps; 3 ns; 6 us;
biggest_corr_time=5 # same as above
N_exp_to_fit=100 # number of exponential functions to be fitted between the samlles and biggest corr time
analyze=1/50 # the proportin of correlation data to be used for fitting, ex. 1/2 uses first half of the data
magnetic_field=2.35 # 5.99 # 8.49 T (values used in SDS paper, J.Chem. Soc.,, Faraday Trans. 1, 1988, 84(12), 4475-4486)
magn_field=magn_field
magnetic_field=magn_field*2*np.pi/gammaH*10**6
nuclei="15N" #nuclei to calculate: 2H-deutherium; 13C - carbon; 15N - nitrogen 



##############3
## CHANGE IN THE CODE 6.4.2022, not going throught the whole content of the folder anymore
###############
take_all_in_folder="number" #"yes"/"no"/"number" analyze all in folder? useful for proteins, if no, fill the following line, if yes fill the folder path
input_corr_file="alphaCF.xvg"

input_prefix="NHrotaCF_" # mostly for peptides, works with take_all_in_folder="no"



## eElab 31.5.22
folder_path= os.getcwd()
folder_path=folder_path + "/correlation_functions/"
directories = folder_path.split("/")
output_name="tst.out"

if "results" in folder_path:
	base  = "/".join(directories[:-2])
	base = base.replace("/results", "")
else:
	base  = "/".join(directories[:-4])


with open(base + '/model_01.pdb', 'r') as file:
	lines = file.readlines()
	words = lines[-3].split()	
	res_nr=words[5]
	print(res_nr)	
sys.path.insert(1, os.path.dirname(base) + '/simulation_scripts/MD_scripts/relaxation_times')
import relaxation_times as rt

author_name="Samuli Ollila"


# In[3]:


import os

os.remove('Ctimes_Coeffs.csv') if os.path.exists('Ctimes_Coeffs.csv') else None

"""Execute the code - this part needs not be modified"""
#rt.initilize_output(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name,author_name)
if take_all_in_folder=="yes":
    for file in os.listdir(sorted(folder_path)):
        with open('Ctimes_Coeffs.txt', 'a') as f:  # Use 'a' (append) mode to add content to the file
            f.write('Res_nr_' + str(res_count) + '\n')
            res_count += 1
        input_corr_file = folder_path+os.fsdecode(file)
        rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name)
elif take_all_in_folder=="number":
#	step_exp=(biggest_corr_time-smallest_corr_time)/N_exp_to_fit
#	Ctimes = 10 ** np.arange(smallest_corr_time, biggest_corr_time, step_exp)
#	Ctimes = Ctimes * 0.001 * 10 ** (-9);
#	Ctimes_to_save=np.zeros([len(Ctimes),residues+1])
#	Ctimes_to_save[:,0]=Ctimes
	with open('relaxation_times.csv', 'w', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(["Residue_nr", "R1_exp", "R1_sim", "R1_diff", "R2_exp", "R2_sim", "R2_diff", "hetNOE_exp", "hetNOE_sim", "hetNOE_diff", "Tau_eff_area"])
	for i in range(1, int(res_nr)+1):
		try:	
			with open('Ctimes_Coeffs.csv', 'a', newline='') as file:
				writer = csv.writer(file)
				writer.writerow(["Res_nr" + str(i) + "_Ctimes_ns", "Res_nr" + str(i) + "_Coeffs"])
			input_corr_file = folder_path+input_prefix+str(i)+".xvg"
			rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name, i)
		except:
			with open('relaxation_times.csv', 'a', newline='') as file:
				writer = csv.writer(file)
				writer.writerow([i, "n", "n", "n", "n", "n", "n", "n", "n", "n", "n"])
else:
    rt.GetRelaxationData(OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,magnetic_field,input_corr_file,nuclei,output_name)


# In[ ]:
