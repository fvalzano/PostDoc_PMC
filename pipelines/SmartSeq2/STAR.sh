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
FASTQ_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/fastq"
OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"
GENOME_DIR=/hpc/pmc_kool/fvalzano/Ref_genome_star/star         # Path to the directory where the STAR genome index is stored
THREADS=8                              # Number of threads to use for STAR alignment
SHORT_CELL_IDS=()
CELL_IDS=()
for FILE in ${FASTQ_DIR}/*/*_R1.fastq.gz; do
    # Extract the base name (e.g., 210FH-P2-A01_R1_trimmed.fq.gz)
    BASENAME=$(basename "$FILE")
    # Remove the '_R1.fastq.gz' suffix to get the ID (e.g., 210FH-P2-A01)
    CELL_ID=$(echo "$BASENAME" | sed 's/_R1.fastq.gz//')
    CELL_IDS+=("$CELL_ID")
done
# Loop through each cell and run STAR alignment

for CELL in "${CELL_IDS[@]}"; do
        # Define input FASTQ files
        FASTQ_R1="${FASTQ_DIR}/*/${CELL}_R1.fastq.gz"
        FASTQ_R2="${FASTQ_DIR}/*/${CELL}_R2.fastq.gz"
    
        # Define output directory for each cell
        CELL_OUTPUT_DIR="${OUTPUT_DIR}/${CELL}"
        mkdir -p ${CELL_OUTPUT_DIR}
    
        # Run STAR alignment in TranscriptomeSAM mode - as we will use subread to count the features
        STAR --runThreadN ${THREADS} \
             --genomeDir ${GENOME_DIR} \
             --readFilesIn ${FASTQ_R1} ${FASTQ_R2} \
             --readFilesCommand zcat \
             --outFileNamePrefix ${CELL_OUTPUT_DIR}/${CELL}_ \
             --outSAMtype BAM Unsorted \
             --quantMode TranscriptomeSAM \
             --outSAMunmapped Within \
             --outSAMattributes Standard

        echo "Alignment completed for ${CELL}"
    done
echo "All alignments are completed!"