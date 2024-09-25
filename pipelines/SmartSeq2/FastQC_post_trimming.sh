#!/bin/bash

#SBATCH --time=10:00:00 # Walltime
#SBATCH --mem=10G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=fastqc_post_trimming

BASE_DIR=/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/trimmed_fastq
OUTPUT_DIR=/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastqc_reports/fastqc_reports_post_trimming
mkdir -p $OUTPUT_DIR

module use --append /hpc/local/Rocky8/pmc_research/etc/modulefiles
module load fastqc/0.11.9

# Find all FASTQ files recursively in the base directory
# Use -name '*.fq.gz' - trim_galore produce files with extension fq.gz
find $BASE_DIR -type f \( -name '*.fq.gz' \) | while read file; do
    echo "Running FastQC on $file..."
    # Run FastQC and save the reports in the output directory
    fastqc "$file" -o $OUTPUT_DIR
done