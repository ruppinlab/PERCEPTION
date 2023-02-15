# PERsonalized single-Cell Expression-based Planning for Treatments In ONcology (PERCEPTION)
We build a precision oncology computational approach capitalizes on recently published matched bulk and single-cell (SC) transcriptome profiles of large-scale cell-line drug screens to build treatment response models from patients' SC tumor transcriptomics. The general objective of this project is to utilize single-cell omics from patients tumor to predict response and resistance. The following figure describe the pipeline of PERCEPTION and its application.
![PERCEPTION pipeline_new](https://user-images.githubusercontent.com/26137763/216864419-a8a35147-ce34-43bd-a360-38d9ac07d4ac.png)

<h2><b>Section 1: How to reproduce the Figures from the PERCEPTION manuscript?</b></h2>
This section will describe how to reproduce the figures described in the following manuscript.
<br><br>
Sinha, S., Vegesna, R., Mukherjee, S., Kammula, A., Dhruba, S.R., Wu, W., Kerr, L., Nair, N., Jones, M., Yosef, N., Stroganov, O., Grishagin, I., Aldape, K., Blakely, C., Jiang, P., Thomas, C., Benes, C., Bivona, T., Schäffer, A. and Ruppin, E. (2023). "<b>PERCEPTION: Predicting patient treatment response and resistance via single-cell transcriptomics of their tumors</b>".
<br><br>
The code to generate the Figures is in R and is presented as a set of R markdown(.Rmd) files to be run in Rstudio. 
The data files are in a mixture of RDS, csv, tsv, and other formats.
Because we want to have version control and quickly do updates, the code is in a GitHub repository.
Because the data files are large and because journals require a frozen snapshot of data, the data files are in a Zenodo repository.
To reproduce the results one must clone the GitHub repository and download the files from the Zenodo repository into a common directory on a computer running Linux or MacOS. We have tested the following instructions on one computer of each type. 


<h3><b>Please perform the following steps to reproduce the Figures from the manuscript:</b></h3>

<b>Step 1:</b> Clone/Download this repository using the command
```
	git clone https://github.com/ruppinlab/SCPO_submission
```
This git clone command should create a directory called SCPO_submission that contains this README.md file and two subdirectories called Tools and Figures.
The Tools subdirectory contains the .Rmd files to be run.
The Figures subdirectory contains the Figures. 

<b>Steps 2:</b> Run the following four unix commands
```
     mkdir Figures_expected 
     cd Figures
     cp * ../Figures_expected
     cd ..
```
These commands copy the figures provided into a new folder called Figures_expected, so that if you run the Figure generating commands
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
<b>Step 6:</b> In each of 
	Step2_Figure1.Rmd  
	Step3_Figure3.Rmd  
	Step5_Figure5_and_Extended_Figure9_10.Rmd  
	Step6_ExtendedFigure4-6.Rmd  
	Step7_ExtendedFigure1.Rmd
     
you need to edit the file once manually to change a line near the top that currently reads: working_directory='/Users/mukherjees11/Documents/SCPO_submission/' Replace the characters in single quotes with the full path to the repository SCPO_submission.
 

<b>Step 7:</b> From the SCPO directory open RStudio 
        You may need to run setwd() to set the working directory.

<b>Step 8:</b> "Step0A" comprises a list of libraries from line 372-399 that are required to run the following scripts. Please install them in your Rstudio environment. [Please refer to this to install a new library: https://r-coder.com/install-r-packages/]

<b>Step 9:</b> Using the File pull-down menu in Rstudio, load and run in succession:
    Tools/Step2_Figure1.Rmd;  Tools/Step3_Figure3.Rmd;  Tools/Step5_Figure5_and_Extended_Figure9_10.Rmd;  Tools/Step6_ExtendedFigure4-6.Rmd;  Tools/Step7_ExtendedFigure1.Rmd.
The filenames are chosen so that  "StepN" represents the order of the results are produced and presented in the manuscript.

<b>Step 10</b> (Optional, open-ended): Using whatever criteria you wish compare the newly generated figures in the Figures subdirectory to the manuscript figures
    in the Figures_expected subdirectory

<b>Notes on the structure of repository:</b>
1. Each script aims (Except Step0A and Step0B) to reproduce at least one figure from either main text or supplementary. 
2. Every script is self-contained.
3. The scripts other than Step0A and Step0B and Step1 can be run in any order.
4. The two scripts starting with "step0" provides the necessary functions to run the rest of the scripts and are accordingly called in the beginning of each scripts.
5. "Step1" is essential to build the PERCEPTION model. By default, this step uses the pan-cancer to build the response model for 44 drugs. However, users can modify this script to change the ‘cancer type’ to build the model for individual cancer types and also can run the PERCEPTION model for specific drugs. However, we found that the pan-cancer model shows better performances than the model built based on individual cancer types.

Scripts to expect from Github repo: A folder name "Tools" comprising following files: 

	step0A_functions_needed.Rmd 
	step0B_functions_needed.Rmd
	Step1_Build_Perception_models.Rmd 
	Step2_Figure1.Rmd 
	Step3_Figure3.Rmd 
	Step5_Figure5_and_Extended_Figure9_10.Rmd
	Step6_ExtendedFigure4-6.Rmd 
	Step7_ExtendedFigure1.Rmd

<b>Files to get from zenodo and place in the Data subdirectory:</b> 

	 carfilzomib_lenalidomide_Bulk_models_list.RDS 
	 carfilzomib_lenalidomide_model.RDS 
	 carfilzomib_lenalidomide_model_bulk_rf.RDS 
	 DepMapv12.RDS 
	 EGFR_WT_signature.csv 
	 EGFR_WT_signature_PDC.csv 
	 FDA_approved_drugs_models.RDS 
	 genes_across_scRNA_datasets_ofInterest.RDS 
	 genesUsed_toBuild.RDS 
	 lung_killing_n_abundance_clone_keyDrugs.RDS 
	 lung_tSNE (1).txt
	 Maynard_Supp2_Demo.xlsx 
	 model_performances.RDS 
	 Predicted_viability_ScreenNishanthApr17.RDS 
	 PRJNA591860.RDS
	 PRJNA591860_patient_level_killing.RDS 
	 PRJNA591860_sample_cell_names.RDS 
	 Responde_models_usedin_PRJNA591860_lung.RDS
	 Response_model_nutlin3.RDS
	 Ruppin120321.xlsx 
	 subClone_id_PRJNA591860.RDS 
	 summary_AUC.tsv 
	 summary_combAUC.tsv 
	 Supp_COhen_IdoAmit_etal.zip 
	 Suppl_table_processed_drugresponse.xlsx 
	 TableS5.csv 
	 TableS5 copy.csv 

<h2><b>Section 2: How to utilize PERCEPTION for a new clinical trial dataset with sc-expression?</b></h2>
To provide guidance to the users on how PERCEPTION could be used for new clinical trial datasets, we have generated a new script named <b>“Running_PERCEPTION_for_new_dataset.Rmd”</b>. This script will provide the fundamental functionality of PERCEPTION in predicting treatment response using a patient's scRNA-seq expression data for particular drugs or combination therapies. 
<br><br>
It should be emphasized that the application of PERCEPTION is not limited to the aforementioned script. Researchers can utilize the PERCEPTION model to explore various aspects of their studies. The scripts (Step2_Figure1.Rmd , Step3_Figure3.Rmd, Step5_Figure5_and_Extended_Figure9_10.Rmd, Step6_ExtendedFigure4-6.Rmd, Step7_ExtendedFigure1.Rmd) described in the Section 1 has already provided examples of this, where we employed PERCEPTION to analyze three different clinical trial datasets for different cancer types and addressed several aspects associated with individual datasets.

<h4><i>The guidance for using "Running_PERCEPTION_for_new_dataset.Rmd" script.</i></h4>
To run the PERCEPTION on the new clinical trial dataset, please use <b>Tools/Running_PERCEPTION_for_new_dataset.Rmd</b>. 
<br><br>
The two main steps are Step1 and Step2 in this script, where attention is needed. In Step1 (lines 32-38), the users need to upload their own files; in Step2 (line 50), the users need to provide the name of their drug of interest for which they want to build the PERCEPTION model and predict the treatment response.
<br><br>
We provided an example of how the script utilizes sample datasets to build the PERCEPTION model for two drugs and subsequently predict treatment responses in a particular type of cancer.
<br>
<h4><i>How to generate the input files for the new dataset?</i></h4>
We have tested this code for a demo dataset; the authors should use this example to format their input files. In Step1 (lines 32-38) of Running_PERCEPTION_for_new_dataset.Rmd, the example input files are provided. The authors should follow these file types to build their own datasets. Basically, three input files are needed as input (i) single-cell gene expression matrix from patients, (ii) cell names in each patient, and (iii) patient demographics, including clinical response data.

<h4><i>Where are the test datasets available?</i></h4>
We used the PRJNA591860 dataset as an example. We are unable to upload these files in GitHub due to restrictions on file size. Please download the PRJNA591860.RDS and PRJNA591860_sample_cell_names.RDS from the Zenodo repository https://doi.org/10.5281/zenodo.6528990 and placed it into the Data directory. If you already downloaded all the files from Zenodo and placed them in the Data directory; then you don't need to download these two files again. PRJNA591860.RDS contains the single-cell gene expression matrix and PRJNA591860_sample_cell_names.RDS contains the cell names in each patient. 
<br><br>
For patient demographics and clinical response data, we have uploaded a sample file named Sample_data_response.xlsx in the Test/Test_data directory. Please download this file in your Data directory.

<h4><i>How to build the model and run the code for specific drugs?</i></h4>
The Step2 of this code (Running_PERCEPTION_for_new_dataset.Rmd) demonstrated how the users could use it to build the response model for specific drugs. We have used the example of two drugs, dabrafenib and erlotinib. The users could select any combination of drugs from the list of 44 drugs to build their own response model. We have described it within the code.

<h4><i>What are the expected outputs of Running_PERCEPTION_for_new_dataset.Rmd?</i></h4>
This script identifies the major cancer cell clusters in the patient's tumor using their sc-expression (transcriptional clone). By computing the mean expression of each transcriptional clone and providing it as an input to the drug response model, this script predicts the drug response of each transcriptional clone separately. Considering that the most resistant clone to the drug will likely get selected by the treatment, this script predicts the overall patient’s response as the predicted response of the most resistant clone.
<br><br>
The figures generated using the sample dataset are provided in the Test/Test_figures directory. 

<br><br>
<b>Contact:</b> Sanju Sinha (sanju@terpmail.umd.edu) 
<br>
Cancer Data Science Laboratory (CDSL), National Cancer Institute (NCI), National Institutes of Health (NIH), Bethesda, Maryland, USA
