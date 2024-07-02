library(readxl)
library(dplyr)
library(stringr)
#Load Overview files containing different IDs
Overview = read_xlsx("/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Overview/General_overview_patient_seqdata.xlsx")
#Delete potential space in the PMCID - this creates troubles during the fetching procedure
Overview=Overview %>% 
  mutate(across(where(is.character), str_trim))
#Specify PMCIDs of interest
IDs = c("806AAS", "222AAS", "745AAS")
IDs = IDs[order(IDs, decreasing = F)]
#Fetch the different IDs for sample type
IDs_RNA_Patient_Biosource= Overview[Overview$PMCID %in% IDs,]$'RNAseq Patient Biosource'
IDs_RNA_Tumoroid_Biomaterial= Overview[Overview$PMCID %in% IDs,]$'RNAseq Tumoroid Biomaterial'
writeLines(IDs_RNA_Patient_Biosource, "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/IDs/RNAseq_Patient_Biosource.txt")
writeLines(IDs_RNA_Tumoroid_Biomaterial, "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/IDs/RNAseq_Tumoroid_Biomaterial.txt")
