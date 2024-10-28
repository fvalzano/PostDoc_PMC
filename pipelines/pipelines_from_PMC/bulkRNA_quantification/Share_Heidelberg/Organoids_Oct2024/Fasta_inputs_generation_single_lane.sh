#!/bin/bash
#SBATCH -o /hpc/pmc_kool/Bulk_RNA/Mieke_organoids/90-1077929049/log.out

# Base directory one level above the directory containing the Fastq files
dir="/hpc/pmc_kool/Bulk_RNA/Share_Heidelberg/Organoids_Oct2024/view-by-pid/"
fq_dir="Fastq"

# Navigate to the base directory
cd ${dir}/${fq_dir}

# Iterate over each subfolder in Fastq
        for subdir in */; do
            cd "$subdir"
            # Get the sample name
            sample=$(for name in `ls`; do echo $name | cut -d _ -f 1; done | uniq)
            # Get the forward and reverse file paths
            forward=$(realpath $(ls | grep _R1.fastq.gz))
            reverse=$(realpath $(ls | grep _R2.fastq.gz))
            # Write the .fasta.inputs file
            printf "001\t${sample}\t${forward}\t${reverse}" > "${sample}.fasta.inputs"
            # Create directories for each samples - change at need in your favourite directory
            mkdir -p "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs/Share_Heidelberg/Oct2024/"
            # Transfer the newly generated fasta in your favourite directory
            scp "${sample}.fasta.inputs" "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs/Share_Heidelberg/Oct2024/"
            # Clean up - we don't want same files all over the places :)
            rm "${sample}.fasta.inputs"
     # Return to the Fastq directory
            cd ..
        done
