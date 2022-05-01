The README.md is part of the repository to replicate the results/Figures of the manuscript entitled
"Predicting patient treatment response and resistance via single-cell transcriptomics of their tumors". 
The general objective of this project is to utilize single-cell omics from patients tumor to predict response and resistance.
In some places, we use the acronym SCPO (single-cell personalized oncology).

The code to generate the Figures is in R and is presented as a set of R markdown(.Rmd) files to be run in Rstudio. 
The data files are in a mixture of RDS, csv, tsv, and other formats.
Because we want to have version control and quickly do updates, the code is in a GitHub repository.
Because the data files are large and because journals require a frozen snapshot of data, the data files are in a Zenodo repository.
To reproduce the results one must clone the GitHub repository and download the files from the Zenodo repository into a common directory on a computer running Linux or MacOS. We have tested the following
instructions on one computer of each type. 

Please perform the following steps to reproduce the main Figures for the SCPO project:
Step 1: Clone/Download this repository using the command 
git clone https://github.com/ruppinlab/SCPO_submission

This git clone command should create a directory called SCPO_submission that contains this README.md file and two subdirectories called Tools and Figures
The Tools subdirectory contains the .RmD files to be run.
The Figures subdirectory contains the Figures. 

Steps 2: Run the following four unix commands
     mkdir Figures_expected 
     cd Figures
     cp * ../Figures_expected
     cd ..
These commands copy the figures provided into a new folde called Figures_expected, so that if you run the Figue generating commands
and thereby replace the Files in the Figures subdirectory, you can subsequently compare the newly generated fiels in Figures
to the expected files in Figures_expected.

Step 3: Run
       mkdir Data
This command creates a Data subdirectory that is parallel to the Figures and Tools subdirectories.
The Data subdirectory will contain the data files, which you will download from Zenodo in the next step.

Step 4: Download the data files from the following Zenodo repository 10.5281/zenodo.6499899
into the Data subdirectory. Due to the limitations of Zenodo, you must download one file at a time. The data files
expected are listed near the bttom of this file. 

Step 5: While in the Data directory,
       run 
       unzip Supp_COhen_IdoAmit_etal.zip
       cd ../Tools

Step 6: In each of Step2_Figure1.Rmd  Step3_Figure3.Rmd  Step5_Figure5_and_Extended_Figure9_10.Rmd  Step6_ExtendedFigure4-6.Rmd  Step7_ExtendedFigure1.Rmd
you need to edit the file once manually to change a line near the top that currently reads:
working_directory='/Users/sinhas8/SCPO_submission/'
Replace the characters in single quotes with the full path to the repository SCPO_submission

Step 7: From the SCPO directory open RStudio 
        You may need to run setwd() to set the working directory.

Step 8: Using the File pull-down menu in Rstudio, load and run in succession:
    Tools/Step2_Figure1.Rmd  Tools/Step3_Figure3.Rmd  Tools/Step5_Figure5_and_Extended_Figure9_10.Rmd  Tools/Step6_ExtendedFigure4-6.Rmd  Tools/Step7_ExtendedFigure1.Rmd
The filenames are chosen so that  "StepN" represents the order of the results are produced and presented in the manuscript.

Step 9 (Optional, open-ended): Using whatever criteria you wish compare the newly generated figires in the Figures subdirectory to the manuscript figures
    in the Figures_expected subdirectory

Notes on the structure of repository:
1. Each script aims (Except Step0) to reproduce a figure from either main text or supplementary. 
2. Every script is self-contained.
3. Scripts could be run in any order necessary.
4. Scripts starting with "step0" provides the necessary functions to run the rest of the scripts and are accordingly called in the beginning of each scripts.

Scripts to expect from Github repo: A folder name "Tools" comprising following files: step0A_functions_needed.Rmd, step0B_functions_needed, Step2_Figure1, Step3_Figure3, Step5_Figure5_and_Extended_Figure9_10, Step6_ExtendedFigure4-6, Step7_ExtendedFigure1

Files to get from zenodo: A single repository name "Data" comprising following files: PRJNA591860_sample_cell_names.RDS, PRJNA591860.RDS, Supp_COhen_IdoAmit_etal.zip, TableS5 copy.csv, Suppl_table_processed_drugresponse.xlsx, Predicted_viability_ScreenNishanthApr17.RDS, summary_combAUC.tsv, summary_AUC.tsv, TableS5.csv, PRJNA591860_patient_level_killing.RDS, Ruppin120321.xlsx, subClone_id_PRJNA591860.RDS, Maynard_Supp2_Demo.xlsx, lung_killing_n_abundance_clone_keyDrugs.RDS, carfilzomib_lenalidomide_model_bulk_rf.RDS, carfilzomib_lenalidomide_Bulk_models_list.RDS, carfilzomib_lenalidomide_model.RDS, genes_across_scRNA_datasets_ofInterest.RDS, genesUsed_toBuild.RDS, model_performances.RDS, EGFR_WT_signature_PDC.csv, EGFR_WT_signature.csv, lung_tSNE.txt, Response_model_nutlin3.RDS, FDA_approved_drugs_models.RDS, DepMapv12.RDS

Contact: Sunju Sinha (sanju@terpmail.umd.edu)
