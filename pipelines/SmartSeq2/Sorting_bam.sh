#!/bin/bash
#SBATCH --time=24:00:00  # Walltime
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mail-type=ALL
#SBATCH -o /hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/log_sort_bam.out
#SBATCH --job-name=SortBAM

#Reason behind this is that counting reads mapping to features it easier for name-sorted files, as read pairs are next to each other
# Load necessary modules
module use --append /hpc/local/Rocky8/pmc_research/etc/modulefiles
module load samtools

# Define the directory containing BAM files
BAM_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/star_alignments"
SORTED_BAM_DIR="/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/sorted_bam"

# Create output directory for sorted BAM files if it doesn't exist
mkdir -p $SORTED_BAM_DIR

# Loop through all BAM files in the STAR output directory
for bam_file in ${BAM_DIR}/*/*.bam; do
    sample_name=$(basename $bam_file _Aligned.out.bam)  # Extract sample name from BAM filename
    echo "Sorting BAM file by name: $sample_name"

    # Sort the BAM file by name and save it to the sorted BAM directory
    samtools sort -n -@ 3 -o ${SORTED_BAM_DIR}/${sample_name}_sorted_by_name.bam $bam_file

    echo "Sorted BAM file saved: ${SORTED_BAM_DIR}/${sample_name}_sorted_by_name.bam"
done

echo "All BAM files sorted by name."
