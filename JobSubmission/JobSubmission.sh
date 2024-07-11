#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=400GB
#SBATCH --job-name=Annotation_SingleR
#SBATCH -o log.out
#SBATCH -e errlog.out
#SBATCH --mail-user=f.valzano@prinsesmaximacentrum.nl 

module use --append /hpc/local/Rocky8/pmc_kool/modulefiles
module load R
Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/R/TME/Scripts_TIGIT/04_Annotation_Dataset.R
