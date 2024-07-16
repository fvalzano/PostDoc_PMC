library(Seurat)
library(readr)
library(SingleR)
scrna_mb = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Entity_Splitting/scrna_harmony_Medulloblastoma.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =50000, replace=F)]
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
genes=scrna_kaesmann@assays[["RNA"]]@meta.features$feature_name
metafeatures = scrna_kaesmann@assays[["RNA"]]@meta.features
metafeatures$ENSG = rownames(metafeatures)
scrna_matrix = scrna_kaesmann@assays[["RNA"]]@counts
rownames(scrna_kaesmann@assays[["RNA"]]@counts) = ifelse(scrna_matrix@Dimnames[[1]] %in% metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
DefaultAssay(scrna_mb) = "RNA"
annotations = SingleR(test=scrna_mb@assays$RNA$counts,
                             ref=scrna_kaesmann@assays$RNA$counts, 
                             labels=scrna_kaesmann$author_cell_type,
                             de.method="wilcox", de.n = 10)
transfer.anno = as.data.frame(annotations$labels, row.names = rownames(annotations))
transfer.anno$`annotations$labels` = as.factor(transfer.anno$`annotations$labels`)
scrna_mb <- AddMetaData(scrna_mb, transfer.anno, col.name = "kaesmann_label")
write_rds(scrna_mb, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/Seurat_subsets/scrna_mb.rds")