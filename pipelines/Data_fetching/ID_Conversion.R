library(readxl)
library(dplyr)
library(stringr)
#Load Overview files containing different IDs
Overview = read_xlsx("/hpc/pmc_kool/fvalzano/Overview/General_overview_patient_seqdata.xlsx")
#Delete potential space in the PMCID - this creates troubles during the fetching procedure
Overview=Overview %>% 
  mutate(across(where(is.character), str_trim))
#-----------IMPORTANT:Specify PMCIDs of interest-----------
IDs = c()
IDs = IDs[order(IDs, decreasing = F)]
#Fetch the different IDs for sample type
IDs_RNA_Patient_Biosource= Overview[Overview$PMCID %in% IDs,]$'RNAseq Patient Biomaterial'
IDs_RNA_Tumoroid_Biomaterial= Overview[Overview$PMCID %in% IDs,]$'RNAseq Tumoroid Biomaterial'
writeLines(IDs_RNA_Patient_Biosource, "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/IDs/RNAseq_Patient_Biomaterial.txt")
writeLines(IDs_RNA_Tumoroid_Biomaterial, "/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/IDs/RNAseq_Tumoroid_Biomaterial.txt")
