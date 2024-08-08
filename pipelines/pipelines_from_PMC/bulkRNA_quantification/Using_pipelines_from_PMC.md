## This code allows to use the bulkRNA analysis pipeline from PMC - This script has been optimized fo experiments belonging to the ITCCP4
    project='ITCCP4'

## Input file preparation:
First, run the Fasta_inputs_generation_single_lane.sh

    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/wdl_scripts/bulkRNA_quantification/Fasta_inputs_generation_single_lane.sh

In case of mixed runs, meaning some reads are single lane and other multi lanes, delete the output of the single lane script for the multiple lanes run.
E.g. you have sequencing runs A, B, C. A, B are single lane and C is multi lanes, running Fasta_inputs_generation_single_lane.sh will output a .fasta also for C,  however, this will be in a wrong format. After deleting the output of Fasta_inputs_generation_single_lane.sh for file C, run 

    bash Fasta_inputs_generation_multi_lanes.sh 
    
for sample C only. This will output a .fasta file in the right format
Then, manually adjust the entries in the multiple lanes script and run it as many times as the samples on multiple lanes.

After that, it is time to retrieve the barcode information, for this, run the following script:
    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/wdl_scripts/bulkRNA_quantification/Barcodes_search.sh

Now that all the input files are ready, we can start using the pipeline:
In order to use the pipelines we will need to modify the .json input file for each pipeline used. 
First we will create json files to use the fastq to ubam conversion. Generate the json files with the following script:

    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/wdl_scripts/bulkRNA_quantification/Fastq_ubam_workflow_fv_inputs_json_generation.sh
IMPORTANT: Don't forget to fill in the information in the newly generated json files (I left comments in the boxes that need to be adjusted)  ---> Try to automate?

## ------FASTQ TO UBAM------
Now, convert the fastq in ubam format using the modified script to loop through several samples (instead of calling the same
function multiple times or save the same script multiple times). This will allow one bash command to submit all the samples you are querying.

    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/wdl_scripts/bulkRNA_quantification/${project}/run_fastq_ubam_workflow_fv_looped.sh

Second, we will create json file to use the ubam to counts conversion. Generate the json file with the following script:

    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/wdl_scripts/bulkRNA_quantification/Rna_fusions_germline_snv_inputs_no_molgenis_fv_generation.sh
IMPORTANT: Don't forget to fill in the information in the newly generated json files (I left comments in the boxes that need to be adjusted)   ---> Try to automate?


## ------UBAM TO COUNTS------
Now, convert the ubam files in count tables using the modified script to loop through several samples (instead of calling the same
function multiple times or save the same script multiple times). This will allow one bash command to submit all the samples you are querying.

    bash /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/pipelines_from_PMC/bulkRNA_quantification/${project}/run_rna_fusion_workflow_fv_looped.sh