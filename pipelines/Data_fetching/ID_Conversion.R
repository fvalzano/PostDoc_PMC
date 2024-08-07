library(readxl)
library(dplyr)
library(stringr)
#Load Overview files containing different IDs
Overview = read_xlsx("/hpc/pmc_kool/fvalzano/Overview/General_overview_patient_seqdata_ForHPC.xlsx")
#Delete potential space in the PMCID - this creates troubles during the fetching procedure
Overview=Overview %>% 
  mutate(across(where(is.character), str_trim))
#-----------IMPORTANT:Specify PMCIDs of interest-----------
IDs = c()
IDs = IDs[order(IDs, decreasing = F)]
#Fetch the different IDs for sample type
IDs_RNA_Patient_Biomaterial= Overview[Overview$PMCID %in% IDs,]$'RNAseq Patient Biomaterial'
IDs_RNA_Tumoroid_Biomaterial= Overview[Overview$PMCID %in% IDs,]$'RNAseq Tumoroid Biomaterial'
writeLines(IDs_RNA_Patient_Biomaterial, "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Patient_Biomaterial.txt")
writeLines(IDs_RNA_Tumoroid_Biomaterial, "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Tumoroid_Biomaterial.txt")
#If there are relapse data, run:
IDs_RNA_Relapse_Biomaterial= Overview[Overview$PMCID %in% IDs,]$'RNAseq Patient Relapse Biomaterial'
writeLines(IDs_RNA_Relapse_Biomaterial, "/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/IDs/RNAseq_Relapse_Biomaterial.txt")

