#!/usr/bin/python3

import pandas as pd
import numpy as np
import os
import glob


SIM_DIR = os.getcwd()
BASE_DIR=os.path.dirname(SIM_DIR)
RESULTS=BASE_DIR + '/results/'


csv_files=glob.glob(RESULTS + 'aSyn/rep_to_exp_data/Accepted_cases/Avg_corr_matrix.csv')[0]

def create_pairwise_dict_from_csv(csv_file_path):
	df = pd.read_csv(csv_file_path, index_col=0)

	pairwise_dict = {}
	headers = df.columns.tolist()  # Assuming row headers are the same as column headers
	print(headers)
	for i in range(len(headers)):
		for j in range(len(headers)):
			if j < i:  # Ensure j is smaller than i
				key = (headers[i], headers[j])
				value = df.iloc[i, j]
				if key in pairwise_dict:
					pairwise_dict[key].append(value)
				else:
					pairwise_dict[key] = [value]
					
	sorted_pairs = sorted(pairwise_dict.items(), key=lambda item: item[1], reverse=True)
#	for pair in sorted_pairs[:2000]:
#		print(pair)
	return pairwise_dict

# Example usage with a CSV file

pairwise_dict = create_pairwise_dict_from_csv(csv_files)
