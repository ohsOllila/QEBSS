Semi-automatic MD production and comparison to NMR relaxation times
(Download and extract Automation_tool_scripts_LUMI or Automation_tool_scripts_MAHTI, depending on the cluster you are going to use) 


The first step is to generate the initial replicas. Go to https://idpconformergenerator.readthedocs.io/en/latest/installation.html or follow the installation steps below: 
 
`cd Automation_tool_scripts_***/Idpconfgenerator_automation`

`git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator`
`cd IDPConformerGenerator` 
`conda env create -f requirements.yml` 
`cd ..`


Copy your fasta file to Automation_tool_scripts_*** directory and generate replicas by running: 

`conda activate idpconfgen` 
`./create_replicas.sh` (Choose the number of your fasta file) 


Copy Unst_prot, MD_parameter_files, simulation_scripts and env.yml to your project scratch in Mahti or Lumi. Using my Mahti paths as example:

`scp -r Unst_prot MD_parameter_files simulation_scripts env.yml malmcajs@mahti.csc.fi:/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare` 


Log into HPC: 
`ssh malmcajs@mahti.csc.fi` 


Go to folder where you added your simulation folders/scripts: 
`cd /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare`

Add your experimental data file to this folder. Make sure to add the first line with the magnetic field strength in MHz. Any missing data should type "n". Name the file Unst_prot_exp_data.txt 

An example of what it should look like can be seen in the file Unst_hydrolase_exp_data.txt or Unst_alphasynuclein_exp_data.txt

Set up the environment here too:

In Mahti:
 
`module purge`
`module load tykky`
`mkdir env`
`conda-containerize new --prefix env env.yml`

In Lumi:

`module purge`
`module load LUMI`
`module load lumi-container-wrapper`
`mkdir env`
`conda-containerize new --prefix env env.yml`


Go to simulation_scripts/run_dir. This is the folder from which you run everything.

In Mahti nothing needs to be changed. The scripts will ask which one of your available project resources you want to use. Choose the number of the correct project. 

❗️In Lumi you need to change the project number in the slurm scripts manually. Go to simulation_scripts/MD_scripts and change the line `#SBATCH --account=project to #SBATCH --account=project_"your project number"` in files md_prep.sh, md.sh and analysis.sh.

Run scripts in order:

`sh run_prep.sh` 
`sh run_batch_md.sh` 
`sh run_analysis.sh` 



❗️Before you run the next script, you must ensure the previous step was finished. E.g., you can check if the system preparation run was successful by going to /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/Unst_prot and running `ls -lh */*/md*tpr`. Run `ls -lh */*/md*tpr | wc –l` to count the files. If this step was succesful you should have 25 md*tpr files when the jobs are finished. To know if the md run is completed you can similairly type `ls -lh */*/md*gro | wc –l`.This script should be rerun multiple times until you have 25 `md*gro` files. The last script collects the results, which can be found in /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/results/Unst_prot.

OBS: The run_batch_md.sh script prevents double submission of identical jobs so you can rerun this script even though some of your jobs are still running. 
