#!/bin/bash

#SBATCH --time=10:00:00 # Walltime
#SBATCH --mem=10G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=SmartSeq2_Seurat_Setup

module use --append /hpc/local/Rocky8/pmc_kool/modulefiles
module load R/4.3.0

Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/Preprocessing_count_matrices.R