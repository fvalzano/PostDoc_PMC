#!/bin/bash
#SBATCH -o /hpc/pmc_kool/Bulk_RNA/Mieke_organoids/90-1077929049/log.out

# Base directory one level above the directory containing the Fastq files
dir="/hpc/pmc_kool/Bulk_RNA/Mieke_organoids/90-1077929049"
fq_dir="00_fastq"

# Navigate to the base directory
cd ${dir}

# Iterate over each subfolder in the 00_fastq directory (each subfolder is a sample ID)
for sample_dir in ${dir}/${fq_dir}/*; do

    # Check if it is a directory
    if [ -d "${sample_dir}" ]; then
        # Extract sample ID from directory name (assuming directory is named after the sample ID)
        sample=$(basename ${sample_dir})
        # Move into the sample directory
        cd ${sample_dir}
        # Find the R1 FASTQ file in the sample directory (assumes files contain '_R1_001.fastq.gz')
        r1_fastq=$(ls | grep "_R1_001.fastq.gz")
        # Check if the R1 FASTQ file exists
        if [ -n "${r1_fastq}" ]; then
            # Extract barcode from the R1 FASTQ file and store it in a temporary file
            zcat $(realpath ${r1_fastq}) | head -1 | rev | cut -d ":" -f 1 | rev > barcodes_${sample}
            # Create the output directory for barcodes if it doesn't already exist
            mkdir -p "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/barcodes/AZ"
            # Move the barcode file to the designated directory
            scp barcodes_${sample} "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/barcodes/AZ"
            # Clean up by removing the temporary barcode file
            rm barcodes_${sample}
        else
            echo "No R1 FASTQ file found in ${sample_dir}"
        fi
        # Return to the parent directory
        cd ..
    else
        echo "${sample_dir} is not a directory. Skipping."
    fi

done

# Return to the original directory after all processing is done
cd ${dir}
