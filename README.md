Semi-automatic MD production and comparison to NMR relaxation times
(Download and extract Automation_tool_scripts_LUMI recipitory) 


The first step is to generate the initial replicas. Go to https://idpconformergenerator.readthedocs.io/en/latest/installation.html or follow the installation steps below: 
 
cd Automation_tool_scripts_*/Idpconfgenerator_automation

git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator 

cd IDPConformerGenerator 

conda env create -f requirements.yml 

cd ..


Copy your fasta file to Automation_tool_scripts_LUMI directory and generate replicas by running: 

conda activate idpconfgen 
./create_replicas.sh (Choose the number of your fasta file) 


Copy Unst_prot, MD_parameter_files, simulation_scripts and env.yml to your project scratch in LUMI. Using my LUMI paths as example:

scp -r Unst_prot MD_parameter_files simulation_scripts env.yml malmcajs@lumi.csc.fi:/scratch/project_462000285/cmcajsa/systems/forcefield_compare 


Log into Lumi: 
ssh malmcajs@lumi.csc.fi 


Go to folder: 
cd /scratch/project_462000285/cmcajsa/systems/forcefield_compare

Add your experimental data to this folder. Make sure to add the first line with the magnetic field strength in MHz. Any missing data should type "n". 
Name the file Unst_prot_exp_data.txt

An example of what it should look like can be seen in the file Unst_alphasynuclein_exp_data.txt

You can also copy relaxation times T1, T2 and hetNOE from https://bmrb.io/ and run `create_exp.py` to automatically create the file. OBS! Copy only the data, skipping the initial text.


Set up the environment:
 
module purge
module load LUMI
module load lumi-container-wrapper
mkdir env
conda-containerize new --prefix env env.yml

You need to manually add your project number to the scripts. Go to simulation_scripts/MD_scripts:

Change line #SBATCH --account=project in files md_prep.sh, md.sh and analysis.sh


Go to simulation_scripts/run_dir and run scripts in order:

sh run_prep.sh 
sh run_batch_md.sh 
sh run_analysis.sh 


Before you run the next script, you must ensure the previous step was finished. E.g., you can check if the system preparation run was successful by going to /scratch/project_462000285/cmcajsa/systems/forcefield_compare/Unst_prot and running ls -lh */*/md*tpr. Run ls -lh */*/md*tpr | wc –l to count the files. The last script collects the results, which can be found in /scratch/project_462000285/cmcajsa/systems/forcefield_compare/results  
