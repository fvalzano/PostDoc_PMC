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

# Define the directory containing the FASTQ files and output reports from FastQC
FASTQ_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/trimmed_fastq"

# STAR genome directory (make sure the genome is pre-indexed using STAR)
GENOME_DIR="/hpc/pmc_kool/fvalzano/Reference_genomes/refdata-gex-GRCh38-2020-A/star/"

# Define output directory for STAR alignments
ALIGN_OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"
mkdir -p $ALIGN_OUTPUT_DIR

# Number of threads for STAR (adjust as needed)
THREADS=8

# Find all FASTQ files recursively in the base directory - remember, trim_galore outputs trimmed fastq files in .fq.gz format
find $FASTQ_DIR -type f -name '*.fq.gz' | while read fastq_file; do
    # Get the base name of the FASTQ file (without directory and extension)
    sample_name=$(basename "$fastq_file" | sed 's/.fq.gz//;s/.fastq//')
    # Create a directory for this sample's output
    sample_output_dir="$ALIGN_OUTPUT_DIR/$sample_name"
    mkdir -p $sample_output_dir
    echo "Running STAR alignment on $fastq_file..."
    # Run STAR alignment
    STAR --genomeDir $GENOME_DIR \
         --readFilesIn "$fastq_file" \
         --runThreadN $THREADS \
         --outFileNamePrefix "$sample_output_dir/$sample_name" \
         --outSAMtype BAM SortedByCoordinate \
         --readFilesCommand zcat # Only if the files are gzipped, remove if not
    echo "Alignment completed for $fastq_file. Output in $sample_output_dir"
done

echo "STAR alignment completed for all files."
