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

#Define the directory containing the aligned BAM files
ALIGN_OUTPUT_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/sorted_bam"

#Define the GTF annotation file (update the path to your GTF file)
GTF_FILE="/hpc/pmc_kool/fvalzano/Ref_genome_star/genes/genes.gtf"

#Define the output file for the count matrix
OUTPUT_COUNTS="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/count_matrices"
mkdir -p $OUTPUT_COUNTS

#Number of threads to use for featureCounts (adjust as needed)
THREADS=8

#Find all BAM files in the alignment output directory
BAM_FILES=$(find $ALIGN_OUTPUT_DIR -type f -name '*.bam')

# Loop through sorted BAM files and run featureCounts for each sample
for bam_file in ${ALIGN_OUTPUT_DIR}/*_sorted_by_name.bam; do
    # Extract sample name from BAM filename
    sample_name=$(basename $bam_file _sorted_by_name.bam)
    echo "Running featureCounts for sample: $sample_name"
    
    # Run featureCounts on the BAM file
    featureCounts -T 8 \
                  -a $GTF_FILE \
                  -o ${OUTPUT_COUNTS}/${sample_name}_featureCounts.txt \
                  -g gene_id \
                  -t gene \
                  -p \
                  -B \
                  -C \
    $bam_file

    echo "featureCounts completed for $sample_name"
done

echo "All samples processed with featureCounts."