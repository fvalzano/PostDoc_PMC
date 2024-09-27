#!/bin/bash
#SBATCH --time=48:00:00 # Walltime
#SBATCH --mem=128G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=STAR

module use --append /hpc/local/CentOS7/pmc_research/etc/modulefiles
module load STAR

# Define the directory containing subdirectories with FASTQ files
BASE_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/trimmed_fastq"

# Define the output directory for STAR alignments
STAR_OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"
mkdir -p $STAR_OUTPUT_DIR

# Path to the STAR genome directory (index generated with STAR)
GENOME_DIR="/hpc/pmc_kool/fvalzano/Ref_genome_star/star/"

# Define the number of threads for STAR
THREADS=8

# Find all unique sample names based on the file naming pattern *_R1.fq.gz and *_R2.fq.gz - remember, trim_galore outputs files in .fq.gz format
for sample in $(find $BASE_DIR -type f -name '*_R1_trimmed.fq.gz' | sed 's/_R1_trimmed.fq.gz//' | sort | uniq); do
    # Define the paths for R1 and R2 fastq files
    R1_FILE="${sample}_R1_trimmed.fq.gz"
    R2_FILE="${sample}_R2_trimmed.fq.gz"

    # Extract sample name for output naming (e.g., sample_R1_trimmed.fq.gz -> sample)
    SAMPLE_NAME=$(basename $sample)

    # Define output prefix for STAR
    STAR_PREFIX="${STAR_OUTPUT_DIR}/${SAMPLE_NAME}_"

    # Run STAR alignment
    echo "Running STAR on sample ${SAMPLE_NAME}..."
    STAR --runThreadN $THREADS \
         --genomeDir $GENOME_DIR \
         --readFilesIn $R1_FILE $R2_FILE \
         --readFilesCommand zcat \
         --outFileNamePrefix $STAR_PREFIX \
         --outSAMtype BAM SortedByCoordinate 
    echo "Alignment completed for ${SAMPLE_NAME}"
done

echo "STAR alignment completed for all samples."