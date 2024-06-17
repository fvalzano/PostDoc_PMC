library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)

###All
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_mb = read_rds(paste0(wd,"Seurat_subsets/Post_Entity_Splitting/scrna_harmony_Medulloblastoma.rds"))
scrna_mb_list = SplitObject(scrna_mb, split.by = "Dataset")
for (i in names(scrna_mb_list)) {
  DefaultAssay(scrna_mb_list[[i]]) = "RNA"
  scrna_mb_list[[i]] = SCTransform(scrna_mb_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg= SelectIntegrationFeatures(scrna_mb_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg, value = TRUE)
subtract<-which(hvg %in% ribo.genes)
hvg_filtered<-hvg[-subtract]
scrna_mb = merge(x = scrna_mb_list[[1]], y= scrna_mb_list[-1], merge.data = TRUE, project = "Ep") 
scrna_mb = RunPCA (scrna_mb, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_filtered)
scrna_mb = RunHarmony(scrna_mb, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_mb = RunUMAP(scrna_mb, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_mb = FindNeighbors(object = scrna_mb, reduction = "harmony", dims = 1:30)
i = seq(0.1, 1, by = 0.1)
scrna_mb = FindClusters(scrna_mb, resolution = i)
write_rds(scrna_mb, paste0(wd, "TME_TIGIT/Seurat_subsets/scrna_mb.rds"))

###Myeloid
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_myeloid = read_rds(paste0(wd,"Seurat_subsets/Post_Annotation/scrna_immune_myeloid.rds"))
Idents(scrna_myeloid) = "Entity"
scrna_myeloid = subset(scrna_myeloid, idents= ("Medulloblastoma"))
scrna_myeloid_list = SplitObject(scrna_myeloid, split.by = "Dataset")
for (i in names(scrna_myeloid_list)) {
  DefaultAssay(scrna_myeloid_list[[i]]) = "RNA"
  scrna_myeloid_list[[i]] = SCTransform(scrna_myeloid_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_myeloid= SelectIntegrationFeatures(scrna_myeloid_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_myeloid, value = TRUE)
subtract<-which(hvg_myeloid %in% ribo.genes)
hvg_myeloid_filtered<-hvg_myeloid[-subtract]
scrna_myeloid = merge(x = scrna_myeloid_list[[1]], y= scrna_myeloid_list[-1], merge.data = TRUE, project = "Ep") 
scrna_myeloid = RunPCA (scrna_myeloid, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_myeloid_filtered)
scrna_myeloid = RunHarmony(scrna_myeloid, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_myeloid = RunUMAP(scrna_myeloid, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_myeloid = FindNeighbors(object = scrna_myeloid, reduction = "harmony", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
scrna_myeloid = FindClusters(scrna_myeloid, resolution = i)
Idents(scrna_myeloid) = "SCT_snn_res.1.2"
scrna_myeloid = RenameIdents(scrna_myeloid, c("0" = "Activated TAM",
                                              "1" = "Immunosuppressive TAM1",
                                              "2" = "DC",
                                              "3" = "Activated TAM",
                                              "4" = "Proinflammatory Microglia1",
                                              "5" = "Monocytes",
                                              "6" = "Immunosuppressive TAM2",
                                              "7" = "prolif. TAM",
                                              "8" = "Immunosuppressive TAM3",
                                              "9" = "Astrocytes1",
                                              "10" = "Astrocytes2",
                                              "11" = "Immunosuppressive TAM4",
                                              "12" = "Astrocytes3",
                                              "13" = "Proinflammatory Microglia2",
                                              "14" = "Astrocytes4",
                                              "15" = "MDSC",
                                              "16" = "Damaged Cells",
                                              "17" = "Astrocytes5",
                                              "18" = "Astrocytes6"))
scrna_myeloid$annotation_fv_v2 = scrna_myeloid@active.ident
Idents(scrna_myeloid) = "SCT_snn_res.1.2"
scrna_myeloid = RenameIdents(scrna_myeloid, c("0" = "Activated TAM",
                                              "1" = "Immunosuppressive TAM",
                                              "2" = "DC",
                                              "3" = "Activated TAM",
                                              "4" = "Proinflammatory Microglia",
                                              "5" = "Monocytes",
                                              "6" = "Immunosuppressive TAM",
                                              "7" = "prolif. TAM",
                                              "8" = "Immunosuppressive TAM",
                                              "9" = "Astrocytes",
                                              "10" = "Astrocytes",
                                              "11" = "Immunosuppressive TAM",
                                              "12" = "Astrocytes",
                                              "13" = "Proinflammatory Microglia",
                                              "14" = "Astrocytes",
                                              "15" = "MDSC",
                                              "16" = "Damaged Cells",
                                              "17" = "Astrocytes",
                                              "18" = "Astrocytes"))
scrna_myeloid$annotation_fv_v1 = scrna_myeloid@active.ident
write_rds(scrna_myeloid, paste0(wd, "TME_TIGIT/Seurat_subsets/scrna_immune_myeloid_mb.rds"))

###Lymphoid
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_lymphoid = read_rds(paste0(wd,"Seurat_subsets/Post_Annotation/scrna_immune_lymphoid.rds"))
Idents(scrna_lymphoid) = "Entity"
scrna_lymphoid = subset(scrna_lymphoid, idents= ("Medulloblastoma"))
scrna_lymphoid_list = SplitObject(scrna_lymphoid, split.by = "Dataset")
scrna_lymphoid_list$SCPCA_MB = NULL
for (i in names(scrna_lymphoid_list)) {
  DefaultAssay(scrna_lymphoid_list[[i]]) = "RNA"
  scrna_lymphoid_list[[i]] = SCTransform(scrna_lymphoid_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_lymphoid= SelectIntegrationFeatures(scrna_lymphoid_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_lymphoid, value = TRUE)
subtract<-which(hvg_lymphoid %in% ribo.genes)
hvg_lymphoid_filtered<-hvg_lymphoid[-subtract]
scrna_lymphoid = merge(x = scrna_lymphoid_list[[1]], y= scrna_lymphoid_list[-1], merge.data = TRUE, project = "Ep") 
scrna_lymphoid = RunPCA (scrna_lymphoid, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_lymphoid_filtered)
scrna_lymphoid = RunHarmony(scrna_lymphoid, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_lymphoid = RunUMAP(scrna_lymphoid, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_lymphoid = FindNeighbors(object = scrna_lymphoid, reduction = "harmony", dims = 1:20)
i = seq(0.1, 2, by = 0.1)
scrna_lymphoid = FindClusters(scrna_lymphoid, resolution = i)
scrna_lymphoid = PrepSCTFindMarkers(scrna_lymphoid)
DEG = FindAllMarkers(scrna_lymphoid, group.by = "SCT_snn_res.1.2", min.pct=0.1,min.diff.pct=0.15, only.pos=T)
Idents(scrna_lymphoid) = "SCT_snn_res.1.2"
scrna_lymphoid = RenameIdents(scrna_lymphoid, c("0" = "Cytotoxic T Cells1",
                                                "1" = "Cytotoxic T Cells2",
                                                "2" = "Naive T Cells",
                                                "3" = "NK Cells",
                                                "4" = "T helper Cells",
                                                "5" = "Intermediate Cells",
                                                "6" = "Intermediate Cells",
                                                "7" = "Cytotoxic T Cells3",
                                                "8" = "Intermediate Cells",
                                                "9" = "B Cells",
                                                "10" = "T reg Cells",
                                                "11" = "prolif. T Cells"))

write_rds(scrna_lymphoid, paste0(wd, "TME_TIGIT/Seurat_subsets/scrna_immune_lymphoid_mb.rds"))

