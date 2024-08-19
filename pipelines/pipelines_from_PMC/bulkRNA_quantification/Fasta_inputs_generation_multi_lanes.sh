#!/bin/bash
#For sequencing runs performed on multiple lanes, it is necessary to first merge the Fastq files and then compute the .fasta.inputs on the merged fastq files
#To merge the fastq files use
#In this specific case the sequencing runs on multiple lanes are the following: ITCC-P4_MB0072/Fastq/pr01-f01-r01, ITCC-P4_MB0646/Fastq/pp01-f01-r01 & ITCC-P4_MB0646/Fastq/tp01-f01-r01

# Navigate to the directory one level above the directory containing the Fastq files
multiple_lane_sample='ITCC-P4_MB0646'
file='Fastq'
subdir='tp01-f01-r01'

cd /hpc/pmc_kool/ITCCP4/PMC_MB/${multiple_lane_sample}/${file}/${subdir}
  # Extract unique sample names and process each sample - sample variable is important for the naming of the .fasta.input naming
            for sample in $(ls | cut -d '-' -f 1-2 | uniq); do
                # Initialize the counter for the line identifiers
                counter=1
                
                # Find all forward and reverse pairs for the sample
                for forward in $(ls | grep "${sample}.*_R1.fastq.gz"); do
                    reverse=$(echo "$forward" | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
                    
                    if [ -f "$reverse" ]; then
                        #Create lane id 
                        sample_and_lane=$(echo "$forward" | cut -d '_' -f 1)
                        # Create the line identifier (e.g., 001, 002, 003, ...)
                        line_id=$(printf "%03d" $counter)
                        # Create the sample identifier with lane information
                        sample_with_lane="${sample}-${lane}"
                        # Get the full paths for the forward and reverse files
                        forward_path=$(realpath "$forward")
                        reverse_path=$(realpath "$reverse")
                        # Write the entry to the .fasta.inputs file
                        printf "${line_id}\t${sample_and_lane}\t${forward_path}\t${reverse_path}\n" >> "${sample}.fasta.inputs"
                        # Increment the counter
                        counter=$((counter + 1))
                    fi
                done
            done
# Create directories for each samples - change at need in your favourite directory
mkdir -p "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs/${multiple_lane_sample}/${subdir}"
# Transfer the newly generated fasta in your favourite directory
scp "${sample}.fasta.inputs" "/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs/${multiple_lane_sample}/${subdir}"
# Clean up - we don't want same files all over the places :)
rm "${sample}.fasta.inputs"

