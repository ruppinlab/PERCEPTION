# PERsonalized single-Cell Expression-based Planning for Treatments In ONcology (PERCEPTION)
We build a precision oncology computational approach capitalizes on recently published matched bulk and single-cell (SC) transcriptome profiles of large-scale cell-line drug screens to build treatment response models from patients' SC tumor transcriptomics. The general objective of this project is to utilize single-cell omics from patients tumor to predict response and resistance. The following figure describe the pipeline of PERCEPTION and its application.
![PERCEPTION pipeline_new](https://user-images.githubusercontent.com/26137763/216801587-558f4f2a-8517-470e-b650-063081dc316b.png)

This README.md file is part of the repository to replicate the results/Figures of the manuscript entitled<br>
"Predicting patient treatment response and resistance via single-cell transcriptomics of their tumors". 

The code to generate the Figures is in R and is presented as a set of R markdown(.Rmd) files to be run in Rstudio. 
The data files are in a mixture of RDS, csv, tsv, and other formats.
Because we want to have version control and quickly do updates, the code is in a GitHub repository.
Because the data files are large and because journals require a frozen snapshot of data, the data files are in a Zenodo repository.
To reproduce the results one must clone the GitHub repository and download the files from the Zenodo repository into a common directory on a computer running Linux or MacOS. We have tested the following
instructions on one computer of each type. 

<h3><b>Please perform the following steps to reproduce the main Figures from the manuscript:</b></h3>
<b>Step 1:</b> Clone/Download this repository using the command 
```
git clone https://github.com/ruppinlab/SCPO_submission
```
This git clone command should create a directory called SCPO_submission that contains this README.md file and two subdirectories called Tools and Figures.
The Tools subdirectory contains the .RmD files to be run.
The Figures subdirectory contains the Figures. 

<b>Steps 2:</b> Run the following four unix commands
```
     mkdir Figures_expected 
     cd Figures
     cp * ../Figures_expected
     cd ..
```
These commands copy the figures provided into a new folde called Figures_expected, so that if you run the Figure generating commands
and thereby replace the Files in the Figures subdirectory, you can subsequently compare the newly generated files in Figures
to the expected files in Figures_expected.

<b>Step 3:</b> Run
```
       mkdir Data
```
This command creates a Data subdirectory that is parallel to the Figures and Tools subdirectories.
The Data subdirectory will contain the data files, which you will download from Zenodo in the next step.

<b>Step 4:</b> Download the data files from the following Zenodo repository https://doi.org/10.5281/zenodo.6528990
into the Data subdirectory. Due to the limitations of Zenodo, you must download one file at a time. The data files
expected are listed near the bottom of this file. 

<b>Step 5:</b> While in the Data directory,
```
       run 
       unzip Supp_COhen_IdoAmit_etal.zip
       cd ../Tools
```
<b>Step 6:</b> In each of Step2_Figure1.Rmd  Step3_Figure3.Rmd  Step5_Figure5_and_Extended_Figure9_10.Rmd  Step6_ExtendedFigure4-6.Rmd  Step7_ExtendedFigure1.Rmd
you need to edit the file once manually to change a line near the top that currently reads:
working_directory='/Users/sinhas8/SCPO_submission/'
Replace the characters in single quotes with the full path to the repository SCPO_submission

<b>Step 7:</b> From the SCPO directory open RStudio 
        You may need to run setwd() to set the working directory.

<b>Step 8:</b> "Step0A" comprises a list of libraries from line 372-399 that are required to run the following scripts. Please install them in your Rstudio environment. [Please refer to this to install a new library: https://r-coder.com/install-r-packages/]

<b>Step 9:</b> Using the File pull-down menu in Rstudio, load and run in succession:
    Tools/Step2_Figure1.Rmd  Tools/Step3_Figure3.Rmd  Tools/Step5_Figure5_and_Extended_Figure9_10.Rmd  Tools/Step6_ExtendedFigure4-6.Rmd  Tools/Step7_ExtendedFigure1.Rmd
The filenames are chosen so that  "StepN" represents the order of the results are produced and presented in the manuscript.

<b>Step 10</b> (Optional, open-ended): Using whatever criteria you wish compare the newly generated figires in the Figures subdirectory to the manuscript figures
    in the Figures_expected subdirectory

<b>Notes on the structure of repository:</b>
1. Each script aims (Except Step0) to reproduce a figure from either main text or supplementary. 
2. Every script is self-contained.
3. Scripts could be run in any order necessary.
4. Scripts starting with "step0" provides the necessary functions to run the rest of the scripts and are accordingly called in the beginning of each scripts.
5. "Step1" is essential to build the PERCEPTION model. By default, this step uses the pan-cancer to build the response model for 44 drugs. However, users can modify this script to change the ‘cancer type’ to build the model for individual cancer types and also can run the PERCEPTION model for specific drugs. However, we found that the pan-cancer model shows better performances than the model built based on individual cancer types.

Scripts to expect from Github repo: A folder name "Tools" comprising following files: step0A_functions_needed.Rmd, step0B_functions_needed, Step2_Figure1, Step3_Figure3, Step5_Figure5_and_Extended_Figure9_10, Step6_ExtendedFigure4-6, Step7_ExtendedFigure1

Files to get from zenodo: A single repository name "Data" comprising following files: PRJNA591860_sample_cell_names.RDS, PRJNA591860.RDS, Supp_COhen_IdoAmit_etal.zip, TableS5 copy.csv, Suppl_table_processed_drugresponse.xlsx, Predicted_viability_ScreenNishanthApr17.RDS, summary_combAUC.tsv, summary_AUC.tsv, TableS5.csv, PRJNA591860_patient_level_killing.RDS, Ruppin120321.xlsx, subClone_id_PRJNA591860.RDS, Maynard_Supp2_Demo.xlsx, lung_killing_n_abundance_clone_keyDrugs.RDS, carfilzomib_lenalidomide_model_bulk_rf.RDS, carfilzomib_lenalidomide_Bulk_models_list.RDS, carfilzomib_lenalidomide_model.RDS, genes_across_scRNA_datasets_ofInterest.RDS, genesUsed_toBuild.RDS, model_performances.RDS, EGFR_WT_signature_PDC.csv, EGFR_WT_signature.csv, lung_tSNE.txt, Response_model_nutlin3.RDS, FDA_approved_drugs_models.RDS, DepMapv12.RDS, Responde_models_usedin_PRJNA591860_lung.RDS

<h2><b>How to utilize PERCEPTION for a new clinical trial dataset with sc-expression?</b></h2>
To run the PERCEPTION on the new clinical trial dataset, please use <b>Tools/Step_N.Rmd</b>. 

This script identifies the major cancer cell clusters in the patient's tumor using their sc-expression (transcriptional clone). By computing the mean expression of each transcriptional clone and providing it as an input to the drug response model, this script predicts the drug response of each transcriptional clone separately. Considering that the most resistant clone to the drug will likely get selected by the treatment, this script predicts the overall patient’s response as the predicted response of the most resistant clone.

The two main steps are Step1 and Step2 in this script, where attention is needed. In Step1 (lines 32-38), the users need to upload their own files; in Step2 (line 50), the users need to provide the name of their drug of interest for which they want to build the PERCEPTION model and predict the treatment response.

<h4><i>How to generate the input files for the new dataset?</i></h4>
We have tested this code for a demo dataset; the authors should use this example to format their input files. In Step1 (lines 32-38), the example input files are provided. The authors should follow these file types to build their own datasets. Basically, three input files are needed as input (i) single-cell gene expression matrix from patients, (ii) cell names in each patient, and (iii) patient demographics, including clinical response data.

<h4><i>How to build the model and run the code for specific drugs?</i></h4>
The Step2 of this code demonstrated how the users could use it to build the response model for specific drugs. We have used the example of two drugs, dabrafenib and erlotinib. The users could select any combination of drugs from the list of 44 drugs to build their own response model. We have described it within the code.
<br><br>
Contact: Sanju Sinha (sanju@terpmail.umd.edu) 
<br>
Cancer Data Science Laboratory (CDSL), National Cancer Institute (NCI), National Institutes of Health (NIH), Bethesda, Maryland, USA
