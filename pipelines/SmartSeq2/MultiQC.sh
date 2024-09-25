#!/bin/bash

#SBATCH --time=10:00:00 # Walltime
#SBATCH --mem=10G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=MultiQC

# Initialize conda for the shell session
source $(conda info --base)/etc/profile.d/conda.sh

#Activate conda environment ss2
conda activate ss2

#Now run python code
python3 /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/SmartSeq2/MultiQC.py
