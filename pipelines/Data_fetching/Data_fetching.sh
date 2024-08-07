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
Date="20240807"
Requester="Julie"
Output_folder="${Date}_${Requester}"
cd /hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests
mkdir -p "$Output_folder"
cd ..
RNAseq_Patient_Biomaterial_input_file="/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Patient_Biomaterial.txt"
# Read each label in the input file
while IFS= read -r label; do
mkdir "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests$Output_folder"
    # Search for files containing the label in their name
wget -m -r -nd -np -P "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests$Output_folder" --user fvalzano --password Brindisi.2021 -e robots=off "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept=${label}*gene_id.exon.counts.txt
done < "$RNAseq_Patient_Biomaterial_input_file"

cd /hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests
RNAseq_Relapse_Biomaterial_input_file="/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Relapse_Biomaterial.txt"
# Read each label in the input file
while IFS= read -r label; do
mkdir "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests/$Output_folder"
    # Search for files containing the label in their name
wget -m -r -nd -np -P "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests$Output_folder" --user fvalzano --password Brindisi.2021 -e robots=off "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept=${label}*gene_id.exon.counts.txt
done < "$RNAseq_Relapse_Biomaterial_input_file"

#cd /hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests
#mkdir -p "$Output_folder"
#cd ..
#RNAseq_Tumoroid_Biomaterial_input_file="/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Tumoroid_Biomaterial.txt"
## Read each label in the input file
#while IFS= read -r label; do
#mkdir "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Requests/$Output_folder"
#    # Search for files containing the label in their name
#wget -m -r -nd -np -P "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Requests/$Output_folder" --user fvalzano --password Brindisi.2021 -e robots=off "https://files.bioinf.prinsesmaximacentrum.nl/shares/PMCLAB2020-142/" --accept=${label}*gene_id.exon.counts.txt
#done < "$RNAseq_Tumoroid_Biomaterial_input_file"

cp "$0" "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/Script_copies/Requests/$Output_folder/$(basename "$0")"