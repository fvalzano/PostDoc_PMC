library(tidyverse)
library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
#Load model weights
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#Load scrna dataset
scrna = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/Seurat_subsets/scrna_mb.rds")
#Rename idents according to tumor cells or not - based on known cell markers
Idents(scrna) = "Dataset"
#SCPCA_MB contributes minorly to the dataset, as the low number of cells creates problem for the SCT in the lymphoid compartment, deleting the dataset at this point from later analysis is the best choice
scrna = subset(scrna, idents=c("MB_David",
                                     "Piyush_MB",
                                     "Riemondy",
                                     "Vladoiu_Medulloblastoma"))


Idents(scrna) = "SCT_snn_res.0.8"

scrna = RenameIdents(scrna, c("0" = "Tumor Cells",
                              "1" = "Tumor Cells",
                              "2" = "Tumor Cells",
                              "3" = "Tumor Cells",
                              "4" = "Tumor Cells",
                              "5" = "Tumor Cells",
                              "6" = "Tumor Cells",
                              "7" = "Tumor Cells",
                              "8" = "Tumor Cells",
                              "9" = "Tumor Cells",
                              "10" = "Tumor Cells",
                              "11" = "Myeloid Cells",
                              "12" = "Tumor Cells",
                              "13" = "Tumor Cells",
                              "14" = "Tumor Cells",
                              "15" = "Tumor Cells",
                              "16" = "Tumor Cells",
                              "17" = "Tumor Cells",
                              "18" = "Tumor Cells",
                              "19" = "Tumor Cells",
                              "20" = "Tumor Cells",
                              "21" = "Tumor Cells",
                              "22" = "Tumor Cells",
                              "23" = "Tumor Cells",
                              "24" = "Tumor Cells",
                              "25" = "Tumor Cells",
                              "26" = "Tumor Cells",
                              "27" = "Tumor Cells",
                              "28" = "Tumor Cells",
                              "29" = "Mesenchymal Cells",
                              "30" = "Tumor Cells",
                              "31" = "Lymphoid Cells",
                              "32" = "Tumor Cells",
                              "33" = "Tumor Cells",
                              "34" = "Tumor Cells",
                              "35" = "Mesenchymal Cells",
                              "36" = "Tumor Cells",
                              "37" = "Tumor Cells",
                              "38" = "Tumor Cells",
                              "39" = "Tumor Cells",
                              "40" = "Tumor Cells",
                              "41" = "Tumor Cells",
                              "42" = "Tumor Cells",
                              "43" = "Tumor Cells",
                              "44" = "Tumor Cells",
                              "45" = "Tumor Cells",
                              "46" = "Tumor Cells",
                              "47" = "Tumor Cells"))
scrna$Major_classes_FV = scrna@active.ident
#Subset scrna object retaining only cells belonging to Tumor Cells and Mesenchymal Cells
scrna_subset = subset(scrna, idents = c("Tumor Cells", "Mesenchymal Cells"))
#Extract cell barcodes with associated annotation from cells belonging to Tumor Cells and Mesenchymal Cells
Meta_non_immune_cells = as.data.frame(scrna_subset$Major_classes_FV)
colnames(Meta_non_immune_cells) = "annotation"

#Load annotated Myeloid and Lymphoid objects with detailed annotation
scrna_lymphoid = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_lymphoid_mb.rds")
scrna_myeloid = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds")
#Extract cell barcodes with associated annotation from metadata slot with detailed cluster information
Meta_myeloid_cells = as.data.frame(scrna_myeloid$annotation_fv_v2)
colnames(Meta_myeloid_cells) = "annotation"
Meta_lymphoid_cells = as.data.frame(scrna_lymphoid$annotation_fv_v2)
colnames(Meta_lymphoid_cells) = "annotation"

#Create master metadata slot via merging the cell barcodes from scrna_subset with the myeloid and lymphoid metadata -> In this way we will create a metadata slot with Tumor cells, Mesenchymal cells and myeloid and lymphoid cell subsets
Meta_master = rbind(Meta_non_immune_cells, Meta_lymphoid_cells, Meta_myeloid_cells)
#Add Meta_master to full scrna object
scrna = AddMetaData(scrna, Meta_master, "annotation_fv_v2")
DimPlot(scrna, group.by = "annotation_fv_v2", label = F)
rm(scrna_myeloid, scrna_lymphoid, scrna_subset)
#Run LIANA interaction
scrna_split = SplitObject(scrna, split.by = "Subgroup")
rm(scrna)
liana_results = list()
liana_results_aggregate = list()
liana_trunc=list()
#Loop liana interaction estimation in dataset split by MB subgroup - NB: Damaged cells have no importance here, remove them
for (i in names(scrna_split)) {
    Idents(scrna_split[[i]]) = "annotation_fv_v2"
    idents = as.character(unique(scrna_split[[i]]$annotation_fv_v2))
    idents = idents[idents != "Damaged Cells"]
    scrna_split[[i]] = subset(scrna_split[[i]], idents = idents)
    liana_results[[i]] = liana_wrap(scrna_split[[i]]) 
    liana_results_aggregate[[i]]= liana_aggregate(liana_results[[i]])
    write.csv2(liana_results_aggregate[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/LIANA/Liana_interaction_", i,".csv"))
}
#Load liana results and filter
liana_files = list.files("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/LIANA")
liana_results_aggregate = list()
for (i in liana_files) {
   liana_results_aggregate[[i]] = read.csv2(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/LIANA/", i))
}
names(liana_results_aggregate) = sub("\\.csv$", "", names(liana_results_aggregate))
#Heatmap visualization of interaction frequency
liana_trunc=list()
labels = list()
for (i in names(liana_results_aggregate)){
    # only keep interactions concordant between methods
    liana_trunc[[i]] = filter(liana_results_aggregate[[i]], aggregate_rank <= 0.05) # note that these pvals are already corrected
p = heat_freq(liana_trunc[[i]], pallette = c("white", "navy"), name = i) 
print(p)
}

    
    

