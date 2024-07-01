#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=Data_fetching
#SBATCH -o log.out
#SBATCH -e errlog.out

#Run conversion of IDs
module use --append /hpc/local/Rocky8/pmc_kool/modulefiles
module load R/4.3.0
Rscript /hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/ID_Conversion.R

#Set up date and requester
Date="20240701"
Requester="Francesco"
Output_folder="${Date}_${Requester}"
cd Requests
mkdir -p "$Output_folder/RNAseq_Patient_Biosource"
cd ..
RNAseq_Patient_Biosource_input_file="IDs/RNAseq_Patient_Biosource.txt"
# Read each label in the input file
while IFS= read -r label; do
mkdir "Requests/$Output_folder/RNAseq_Patient_Biosource/$label"
    # Search for files containing the label in their name
wget -m -r -nd -np -P "Requests/$Output_folder/RNAseq_Patient_Biosource/$label" --user fvalzano --password Brindisi.2021 -e robots=off "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept=${label}*gene_id.exon.counts.txt
done < "$RNAseq_Patient_Biosource_input_file"

cd Requests
mkdir -p "$Output_folder/RNAseq_Tumoroid_Biomaterial"
cd ..
RNAseq_Tumoroid_Biomaterial_input_file="IDs/RNAseq_Tumoroid_Biomaterial.txt"
# Read each label in the input file
while IFS= read -r label; do
mkdir "Requests/$Output_folder/RNAseq_Tumoroid_Biomaterial/$label"
    # Search for files containing the label in their name
wget -m -r -nd -np -P "Requests/$Output_folder/RNAseq_Tumoroid_Biomaterial/$label" --user fvalzano --password Brindisi.2021 -e robots=off "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept=${label}*gene_id.exon.counts.txt
done < "$RNAseq_Tumoroid_Biomaterial_input_file"