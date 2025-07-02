**Semi-automatic Workflow for the Quality Evaluation-Based Simulation Selection (QEBSS)**

QEBSS is a protocol that generates and evaluates a diverse set of molecular dynamics simulations based on various starting structures and force fields. The simulations are evaluated against NMR relaxation times (R1, R2) and heteronuclear NOE (hetNOE). This helps identify the best simulations for intrinsically disordered proteins (IDPs) or partially disordered proteins. The analysis includes contact maps, distance maps, and backbone correlation maps, which could help understand the dynamic nature of different IDPs and their roles in biological processes.

This is automized version of the QEBSS protocol to find ensembles in best agreement with NMR spin relaxation data using MD simulations for instrinsically disordered proteins (IDPs). QEBSS is originally introduced in
> Sandelin, A., Nencini, R., Yasar, E., Fudo, S., Stratoulias, V., Kajander, T., & Ollila, S. (2025). **QEBSS: Quality Evaluation Based Simulation Selection for analysis of conformational ensembles and dynamics of multidomain proteins**. *ChemRxiv* [https://doi.org/10.26434/chemrxiv-2024-h3pmt-v2](https://doi.org/10.26434/chemrxiv-2024-h3pmt-v2)

The automized version included in this repository is implemented in
> Malm, C., Girych, M., & Ollila, O. H. S. (2025). **Conformational Rigidity Classification in Intrinsically Disordered Proteins via Integrated NMR and MD Simulations**. *ChemRxiv*, [https://doi.org/10.26434/chemrxiv-2025-m5m0p](https://doi.org/10.26434/chemrxiv-2025-m5m0p)

For running these simulations you need access to a supercomputer like Mahti or Lumi hosted by [CSC – IT Center for Science](https://csc.fi/en/). Depending on the platform that you are using, use the appropriate folder named Automation_tool_scripts_*. 


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


Copy Unst_prot, MD_parameter_files, simulation_scripts and env.yml to your project scratch in MAHTI/LUMI.

Go to the folder where you copied all the folders and add your experimental data there. Make sure to add the first line with the magnetic field strength in MHz. Any missing data should type "n". 
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


OBS! Before you run the next script, you must ensure the previous step was finished. 
