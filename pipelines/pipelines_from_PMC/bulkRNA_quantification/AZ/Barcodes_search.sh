#!/bin/bash

# Navigate to the directory one level above the directory containing the Fastq files
cd /hpc/pmc_kool/ITCCP4/PMC_MB
for dir in */; do
    # Enter in the Fastq subdirectory of each subfolder
    if [ -d "$dir/Fastq" ]; then
        cd "$dir/Fastq"

# Iterate over each subfolder in Fastq
        for subdir in */; do
            cd "$subdir"
            sample=$(for name in `ls`; do echo $name | cut -d - -f 1-2; done | uniq)  #this extracts the sample identifier
            #Move the barcode file to the input folder 
            zcat $(realpath $(ls | grep _R1.fastq.gz)) | head -1 | rev | cut -d ":" -f 1 | rev > barcodes_${sample}
            # Create directories for each samples - change at need in your favourite directory
            mkdir -p "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/barcodes/${dir}${subdir}"
            # Transfer the newly generated fasta in your favourite directory
            scp barcodes_${sample} "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/barcodes/${dir}${subdir}"
            # Clean up - we don't want same files all over the places :)
            rm barcodes_${sample}
            # Return to the Fastq directory
            cd ..
        done
# Return to the main directory
        cd ../..
    fi
done


