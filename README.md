Repository to replicate the results/Figures of "Predicting patient treatment response and resistance via single-cell transcriptomics of their tumors" project. Utilizing single-cell omics from patients tumor to predict response and resistance.

Please perform the following steps to reproduce results of SCPO project:
Step 1: Clone/Download this repository.

Step 2: Download the needed data to execute these scripts from the following Zenodo repository

Step 3: Move the folder "Data" to "SCPO_submission" folder. Note to unzip following files: "Supp_COhen_IdoAmit_etal.zip", 'cancer_approved_drugs.zip', 'PRJNA591860.zip'

Step 4: Open RStudio and set the "working_directory" variable (top of each script) to the local address where "SCPO_submission" repository is cloned.

Filename Annotation: "StepX" represents the order of the results are produced and presented in the papers. 

Notes on the structure of repository:
1. Each script aims (Except Step0) to reproduce a figure from either main text or supplementary. 
2. Every script is self-contained.
3. Scripts could be run in any order necessary.
4. Scripts starting with "Step0" provides the necessary functions to run the rest of the scripts and are accordingly called in the beginning of each scripts.



