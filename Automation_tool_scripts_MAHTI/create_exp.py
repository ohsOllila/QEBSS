#!/usr/bin/python3

res_nr=72

R1_list=["n"] * res_nr
R2_list=["n"] * res_nr
NOE_list= ["n"] * res_nr



with open("T1_data.txt", 'r') as file:
	i=0
	lines = file.readlines()
	while i < len(lines) and len(lines[i].split()) > 1:
		parts = lines[i].split()
		R1_list[int(parts[4]) - 1] = 1 /float(parts[10])
		i+=1

with open("T2_data.txt", 'r') as file:
	i=0
	lines = file.readlines()
	while i < len(lines) and len(lines[i].split()) > 1:
		parts = lines[i].split()
		R2_list[int(parts[4]) - 1] = 1 / float(parts[10])
		i+=1

with open("NOE_data.txt", 'r') as file:
	i=0
	lines = file.readlines()
	while i < len(lines) and len(lines[i].split()) > 1:
		parts = lines[i].split()
		NOE_list[int(parts[4]) - 1] = float(parts[19])
		i+=1	

with open("Unst_prot_exp_data.txt", 'w') as file:
	file.write("R1 R2 NOE\n")
	for i in range(res_nr):
		res = i + 1
		R1 = R1_list[i]
		R2 = R2_list[i]
		NOE = NOE_list[i]
		file.write(f"{res}\t{R1}\t{R2}\t{NOE}\n")
