#!/bin/bash
#SBATCH --time=48:00:00 # Walltime
#SBATCH --mem=256G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=STAR

module use --append /hpc/local/CentOS7/pmc_research/etc/modulefiles
module load STAR

# Define base directory and output directories
BASE_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastq"
MERGED_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/merged_fq"
STAR_OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"

# Create output directories if they don't exist
mkdir -p $MERGED_DIR
mkdir -p $STAR_OUTPUT_DIR
# Path to the STAR genome directory (index generated with STAR)
GENOME_DIR="/hpc/pmc_kool/fvalzano/Ref_genome_star/star/"

# Define the number of threads for STAR
THREADS=8

# Loop through unique sample names
# Assuming file names follow the pattern: <sample_name>_A01_R1.fastq.gz, <sample_name>_A01_R2.fastq.gz, etc.
for file in $(ls ${BASE_DIR}); do
    cd $file
    for sample in $(ls ${BASE_DIR}/${file} | grep -o '^[^-]*-[^-]*' | sort -u); do
        echo "Processing sample: $sample"

        # Merge R1 files across lanes
        cat $(ls $BASE_DIR/$file/${sample}*_R1.fastq.gz) > $MERGED_DIR/${sample}_merged_R1.fastq.gz
        # Merge R2 files across lanes
        cat $(ls $BASE_DIR/$file/${sample}*_R2.fastq.gz) > $MERGED_DIR/${sample}_merged_R2.fastq.gz

        # Run STAR alignment
        echo "Running STAR alignment for $sample..."
        STAR --genomeDir $GENOME_DIR \
             --readFilesIn $MERGED_DIR/${sample}_merged_R1.fastq.gz $MERGED_DIR/${sample}_merged_R2.fastq.gz \
             --readFilesCommand zcat \
             --runThreadN $THREADS \
             --outFileNamePrefix $STAR_OUTPUT_DIR/${sample}_ \
             --outSAMtype BAM SortedByCoordinate

        echo "STAR alignment completed for $sample."
    cd ..
    done
done
echo "All samples processed."