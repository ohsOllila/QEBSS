**Semi-automatic Workflow for Quality Evaluation-Based Simulation Selection (QEBSS)**

QEBSS is a tool that generates and evaluates a diverse set of molecular dynamics simulations based on various starting structures and force fields. The simulations are compared using NMR relaxation times (R1, R2) and heteronuclear NOE (hetNOE). This helps identify the best simulations for intrinsically disordered proteins (IDPs) or partially disordered proteins.

The analysis includes contact maps, distance maps, and rigidity, which help find parts of the protein that show more coordinated movements, such as hairpin-like structures.


For running these simulations you need access to a supercomputer like Mahti or Lumi. Depending on the platform you are gonna use download and extract the right recipitory named Automation_tool_scripts_*. 


The first step is to generate the initial conformers. Go to https://idpconformergenerator.readthedocs.io/en/latest/installation.html or follow the installation steps below: 
 
`cd Automation_tool_scripts_**/Idpconfgenerator_automation`
`git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator` 
`cd IDPConformerGenerator` 
`conda env create -f requirements.yml` 
`cd ..`


Copy your fasta file to Automation_tool_scripts_*/Idpconfgenerator_automation directory and generate replicas by running: 

`conda activate idpconfgen` 
`./create_replicas.sh` (Choose the number of your fasta file) 

This step will generate five initial structures that you can find in the folder Automation_tool_scripts_*/Unst_prot


Copy Unst_prot, MD_parameter_files, simulation_scripts and env.yml to your project scratch in MAHTI/LUMI. Using my MAHTI paths as example:

`scp -r Unst_prot MD_parameter_files simulation_scripts env.yml malmcajs@mahti.csc.fi:/scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare` 


Log into Mahti: 
`ssh malmcajs@mahti.csc.fi` 

Go to the folder where you copied all the files: 
`cd /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare`

Add your experimental data to this folder. Make sure to add the first line with the magnetic field strength in MHz. Any missing data should type "n". 
Name the file Unst_prot_exp_data.txt. 

An example of what it should look like can be seen in the file Unst_alphasynuclein_exp_data.txt

You can also copy relaxation times T1, T2 and hetNOE from https://bmrb.io/ and run `create_exp.py` to automatically create the file. OBS! Copy only the data, skipping the initial text.


Set up the environment:
 
`module purge`
`module load tykky`
`mkdir env`
`conda-containerize new --prefix env env.yml`


You need to manually add your project number to the scripts in LUMI, in MAHTI you can select the project number when running the scripts. Go to simulation_scripts/MD_scripts:

Change line #SBATCH --account=project in files md_prep.sh, md.sh and analysis.sh


Go to simulation_scripts/run_dir and run scripts in order:

`sh run_prep.sh` 
`sh run_batch_md.sh` 
`sh run_analysis.sh` 


Before you run the next script, you must ensure the previous step was finished. E.g., you can check if the system preparation run was successful by going to /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/Unst_prot and running ls -lh */*/md*tpr. Run ls -lh */*/md*tpr | wc –l to count the files. The last script collects the results, which can be found in /scratch/project_2003809/cmcajsa/MD-stabilization/structures/forcefield_compare/results
