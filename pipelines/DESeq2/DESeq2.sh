#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=64GB
#SBATCH --job-name=DESeq2
#SBATCH -o log.out
#SBATCH -e errlog.out

#Set up date and requester
Date="20240701"
Requester="Francesco"
Output_folder="${Date}_${Requester}"
cd Requests
mkdir -p "$Output_folder"
cd ..

#Run DESeq2
module use --append /hpc/local/Rocky8/pmc_kool/modulefiles
module load R/4.3.0
Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2/DESeq2.R

cp "$0" "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2/Script_copies/Requests/$Output_folder/$(basename "$0")"