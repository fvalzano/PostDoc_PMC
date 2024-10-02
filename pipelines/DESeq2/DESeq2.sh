#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=64GB
#SBATCH --job-name=DESeq2
#SBATCH -o log.out
#SBATCH -e errlog.out

#Set up date and requester
Date="20240807"
Requester="Julie"
Output_folder="${Date}_${Requester}"
cd /hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests
mkdir -p "/hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Requests/$Output_folder"
cd ..

#Run DESeq2
module use --append /hpc/local/Rocky8/pmc_kool/modulefiles
module load R/4.3.0
Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2/DESeq2.R

cd /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2
mkdir -p /hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Script_copies/Requests/$Output_folder
cp DESeq2.R "/hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Script_copies/Requests/$Output_folder/DESeq2.R"