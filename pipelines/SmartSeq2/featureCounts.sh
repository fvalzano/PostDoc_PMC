#!/bin/bash
#SBATCH --time=48:00:00 # Walltime
#SBATCH --mem=128G
#SBATCH --nodes=1          # Use 1 Node     (Unless code is multi-node parallelized)
#SBATCH --ntasks=1        
#SBATCH --cpus-per-task=3 # number of threads we want to run on
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log.out
#SBATCH --job-name=subread

module use --append /hpc/local/CentOS7/pmc_research/etc/modulefiles
module load subread

# Define the directory containing the aligned BAM files
ALIGN_OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"

# Define the GTF annotation file (update the path to your GTF file)
GTF_FILE="/hpc/pmc_kool/fvalzano/Reference_genomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

# Define the output file for the count matrix
OUTPUT_COUNTS="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/count_matrices/"
mkdir -p $OUTPUT_COUNTS


# Number of threads to use for featureCounts (adjust as needed)
THREADS=8

# Find all BAM files in the alignment output directory
BAM_FILES=$(find $ALIGN_OUTPUT_DIR -type f -name '*.bam')

# Run featureCounts on all BAM files to generate the gene count matrix
echo "Running featureCounts on BAM files..."

featureCounts -T $THREADS \
    -a $GTF_FILE \
    -o $OUTPUT_COUNTS \
    -g gene_id \
    -t exon \
    -s 2 \ # Change to 0 if your data is unstranded, or 1 for stranded
    $BAM_FILES

echo "featureCounts completed. Count matrix saved to $OUTPUT_COUNTS."
