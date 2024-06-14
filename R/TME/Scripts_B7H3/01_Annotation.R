library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)

wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_ep = read_rds(paste0(wd,"Seurat_subsets/Post_Entity_Splitting/scrna_harmony_Ependymoma.rds"))
