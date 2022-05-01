Repository to replicate the results/Figures of "Predicting patient treatment response and resistance via single-cell transcriptomics of their tumors" project. Utilizing single-cell omics from patients tumor to predict response and resistance.

Please perform the following steps to reproduce results of SCPO project:
Step 1: Clone/Download this repository.

Step 2: Download the needed data to execute these scripts from the following Zenodo repository

Step 3: Move the folder "Data" to "SCPO_submission" folder. Note to unzip following files: "Supp_COhen_IdoAmit_etal.zip".

Step 4: Open RStudio and set the "working_directory" variable (top of each script) to the local address where "SCPO_submission" repository is cloned.

Filename Annotation: "StepX" represents the order of the results are produced and presented in the papers. 

Notes on the structure of repository:
1. Each script aims (Except Step0) to reproduce a figure from either main text or supplementary. 
2. Every script is self-contained.
3. Scripts could be run in any order necessary.
4. Scripts starting with "Step0" provides the necessary functions to run the rest of the scripts and are accordingly called in the beginning of each scripts.

Scripts to expect from Github repo: A folder name "Tools" comprising following files: step0A_functions_needed.Rmd, Step0B_functions_needed, Step2_Figure1, Step3_Figure3, Step5_Figure5_and_Extended_Figure9_10, Step6_ExtendedFigure4-6, Step7_ExtendedFigure1

Files to get from zenodo: A single repository name "Data" comprising following files: PRJNA591860_sample_cell_names.RDS, PRJNA591860.RDS, Supp_COhen_IdoAmit_etal.zip, TableS5 copy.csv, Suppl_table_processed_drugresponse.xlsx, Predicted_viability_ScreenNishanthApr17.RDS, summary_combAUC.tsv, summary_AUC.tsv, TableS5.csv, PRJNA591860_patient_level_killing.RDS, Ruppin120321.xlsx, subClone_id_PRJNA591860.RDS, Maynard_Supp2_Demo.xlsx, lung_killing_n_abundance_clone_keyDrugs.RDS, carfilzomib_lenalidomide_model_bulk_rf.RDS, carfilzomib_lenalidomide_Bulk_models_list.RDS, carfilzomib_lenalidomide_model.RDS, genes_across_scRNA_datasets_ofInterest.RDS, genesUsed_toBuild.RDS, model_performances.RDS, EGFR_WT_signature_PDC.csv, EGFR_WT_signature.csv, lung_tSNE.txt, Response_model_nutlin3.RDS, FDA_approved_drugs_models.RDS, DepMapv12.RDS


