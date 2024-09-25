#!/bin/bash

#SBATCH --time=10:00:00 # Walltime
#SBATCH --mem=10G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=Filbin_fastqc

BASE_DIR=/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/
OUTPUT_DIR=/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastqc_reports
mkdir -p $OUTPUT_DIR

module use --append /hpc/local/Rocky8/pmc_research/etc/modulefiles
module load fastqc/0.11.9
# Find all FASTQ files recursively in the base directory
# Use -name '*.fastq.gz' 
find $BASE_DIR -type f \( -name '*.fastq.gz' \) | while read file; do
    echo "Running FastQC on $file..."
    # Run FastQC and save the reports in the output directory
    fastqc "$file" -o $OUTPUT_DIR
done