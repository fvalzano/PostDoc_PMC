#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --mem=64GB
#SBATCH --job-name=Data_fetching
#SBATCH -o log.out
#SBATCH -e errlog.out

#Run conversion of IDs
Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/ID_Conversion.R

#Set up date and requester
Date="20240701"
Requester="Francesco"
Output_folder="${Date}_${Requester}"
cd Requests
mkdir -p "$Output_folder"

RNAseq_Patient_Biosource_input_file="IDs/RNAseq_Patient_Biosource.txt"
# Read each label in the input file
while IFS= read -r label; do
mkdir "Requests/$Output_folder/$label"
    # Search for files containing the label in their name
wget -P "Requests/$Output_folder/$label" --user fvalzano --password Brindisi.2021 "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept="PMLBM000IZI.gene_id.exon.counts.txt"
done < "$RNAseq_Patient_Biosource_input_file"