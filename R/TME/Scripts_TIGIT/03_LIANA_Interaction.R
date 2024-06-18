library(tidyverse)
library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(cowplot)
#Load model weights
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#Load scrna dataset
wd = paste0(getwd(), "/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/")
scrna = readRDS(paste0(wd, "Seurat_subsets/scrna_mb.rds"))
#Rename idents according to tumor cells or not - based on known cell markers
Idents(scrna) = "SCT_snn_res.0.7"
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
                              "11" = "Tumor Cells",
                              "12" = "Myeloid Cells",
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
                              "26" = "Mesenchymal Cells",
                              "27" = "Tumor Cells",
                              "28" = "Tumor Cells",
                              "29" = "Tumor Cells",
                              "30" = "Lymphoid Cells",
                              "31" = "Tumor Cells",
                              "32" = "Tumor Cells",
                              "33" = "Mesenchymal Cells",
                              "34" = "Tumor Cells",
                              "35" = "Tumor Cells",
                              "36" = "Tumor Cells",
                              "37" = "Tumor Cells",
                              "38" = "Tumor Cells",
                              "39" = "Tumor Cells",
                              "40" = "Tumor Cells",
                              "41" = "Tumor Cells",
                              "42" = "Tumor Cells",
                              "43" = "Tumor Cells",
                              "44" = "Tumor Cells",
                              "45" = "Tumor Cells"))
scrna$Major_classes_FV = scrna@active.ident
#Subset scrna object retaining only cells belonging to Tumor Cells and Mesenchymal Cells
scrna_subset = subset(scrna, idents = c("Tumor Cells", "Mesenchymal Cells"))
#Extract cell barcodes with associated annotation from cells belonging to Tumor Cells and Mesenchymal Cells
Meta_non_immune_cells = as.data.frame(scrna_subset$Major_classes_FV)
colnames(Meta_non_immune_cells) = "annotation"

#Load annotated Myeloid and Lymphoid objects with detailed annotation
scrna_lymphoid = readRDS(paste0(wd, "Seurat_subsets/Post_Annotation/scrna_immune_lymphoid_mb.rds"))
scrna_myeloid = readRDS(paste0(wd, "Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds"))
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
