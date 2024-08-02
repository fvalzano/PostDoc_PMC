#!/bin/bash

# Navigate to the PMC_MB directory
cd /hpc/pmc_kool/ITCCP4/PMC_MB
for dir in */; do
    # Enter the Fastq subdirectory of each subfolder
    if [ -d "$dir/Fastq" ]; then
        cd "$dir/Fastq"

# Iterate over each subfolder in Fastq
        for subdir in */; do
            cd "$subdir"
            sample=$(for name in `ls`; do echo $name | cut -d - -f 1-2; done | uniq)  #this extracts the sample identifier
            
            
            #Move the barcode file to the input folder - PICKUP FROM HERE
            zcat $(realpath $(ls | grep _R1.fastq.gz)) | head -1 | rev | cut -d ":" -f 1 | rev > barcodes_${sample}.txt
            scp barcodes_${sample}.txt /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/barcodes/
            rm barcodes_${sample}.txt
            # Return to the Fastq directory
            cd ..
        done
# Return to the main directory
        cd ../..
    fi
done


