library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)
library(enrichplot)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(scran)
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
scrna_mb = readRDS(paste0(wd, "TME_TIGIT/Seurat_subsets/scrna_mb.rds"))
Idents(scrna_mb) = "SCT_snn_res.0.8"
scrna_myeloid = subset(scrna_mb, idents="11")
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
write_rds(scrna_myeloid, paste0(wd, "TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds"))

#Load reference datasets and run SingleR autoannotation
reference = list.files(paste0(wd,"Annotation_reference/Antunes_GBM/RDS_file"))
seurat_objects=list()
for (i in reference) {
  seurat_objects[[i]] = read_rds(paste0(wd,"Annotation_reference/Antunes_GBM/RDS_file/", i))
  Idents(seurat_objects[[i]]) = "cluster"
}
Ref1_Myeloid = subset(seurat_objects[[1]], idents = c("TAM 2", "Monocytes", "TAM 1", "prol. TAM", "DC"))
Ref2_Myeloid = subset(seurat_objects[[2]], idents = c("DC 1", "DC 2", "DC 3", "DC 4", "Monocytes", "prol. TAM", "TAM 1", "TAM 2"))
Ref3_Myeloid = subset(seurat_objects[[3]], idents = c("Hypoxic Mg-TAM", "IFN Mg-TAM", "Mg-TAM", "Phago/Lipid Mg-TAM"))
Ref1_Lymphoid = subset(seurat_objects[[1]], idents = c("B cells", "NK cells", "T cells"))
Ref2_Lymphoid = subset(seurat_objects[[2]], idents = c("B cells", "NK cells", "Plasma B", "Regulatory T cells", "T cells"))

#Get first level of annotation through SingleR's autoannotation 
scrna_myeloid = read_rds(paste0(wd, "TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds"))
annotations_top = SingleR(test = scrna_myeloid@assays$RNA$counts,
                          ref = list(Ref1_Myeloid@assays$RNA$counts,
                                     Ref2_Myeloid@assays$RNA$counts,
                                     Ref3_Myeloid@assays$RNA$counts),
                          labels = list(Ref1_Myeloid$cluster,
                                        Ref2_Myeloid$cluster,
                                        Ref3_Myeloid$cluster),
                          de.method="wilcox")
transfer.anno = as.data.frame(annotations_top$labels, row.names = rownames(annotations_top))
transfer.anno$`annotations_top$labels` = as.factor(transfer.anno$`annotations_top$labels`)
scrna_myeloid <- AddMetaData(scrna_myeloid, transfer.anno, col.name = "Antunes_annotations_top")
scrna_myeloid = PrepSCTFindMarkers(scrna_myeloid)
Idents(scrna_myeloid) = "SCT_snn_res.1"
DEG_v1 = FindAllMarkers(scrna_myeloid, 
                               only.pos = T, 
                               min.pct = 0.15,
                               min.diff.pct=0.15)
DEG_v1 = DEG_v1[DEG_v1$p_val_adj <= 0.05,]
write.csv2(DEG_v1, paste0(wd, "TME_TIGIT/DEG/Myeloid/DEG_v1.csv"))
#Setup Enrichment analysis 
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
#Split master DEG analysis per cluster and run enrichment analysis
for(i in unique(DEG_v1$cluster)){
  DEG_list[[i]] = DEG_v1[DEG_v1$cluster ==i,]
  DEG_list[[i]] = DEG_list[[i]][order(DEG_list[[i]]$avg_log2FC, decreasing = T),] 
  DEG_list_IDs[[i]] = DEG_list[[i]]$gene
  DEG_list[[i]] = DEG_list[[i]]$avg_log2FC
  DEG_vector[[i]] = unlist(DEG_list[[i]])
  names(DEG_vector[[i]]) = DEG_list_IDs[[i]]
  DEG_vector[[i]] = na.omit(DEG_vector[[i]])
  gse[[i]] <- enrichGO(gene=names(DEG_vector[[i]]),
             ont ="BP", 
             keyType = "SYMBOL", 
             pvalueCutoff = 0.05, 
             minGSSize = 5,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
  write.csv2(gse[[i]], paste0(wd, "TME_TIGIT/Enrichment/Myeloid/v1/Cluster_", i, ".csv"))
}
Idents(scrna_myeloid) = "SCT_snn_res.1"
scrna_myeloid = RenameIdents(scrna_myeloid, c("0" = "Activated TAM",
                                              "1" = "Immunosuppressive TAM",
                                              "2" = "Activated TAM",
                                              "3" = "Astrocytes",
                                              "4" = "Proinflammatory Microglia1",
                                              "5" = "Monocytes",
                                              "6" = "Immunosuppressive TAM",
                                              "7" = "DC",
                                              "8" = "Proinflammatory Microglia2",
                                              "9" = "Hypoxic TAM",
                                              "10" = "Astrocytes",
                                              "11" = "Immunosuppressive TAM",
                                              "12" = "Damaged Cells",
                                              "13" = "Astrocytes",
                                              "14" = "prolif. TAM",
                                              "15" = "DC",
                                              "16" = "Astrocytes",
                                              "17" = "MDSC",
                                              "18" = "Astrocytes",
                                              "19" = "Astrocytes"))
scrna_myeloid$annotation_fv_v2 = scrna_myeloid@active.ident
write_rds(scrna_myeloid, paste0(wd, "TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds"))



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
scrna_lymphoid$annotation_fv_v2 = scrna_lymphoid@active.ident                                                                                   
write_rds(scrna_lymphoid, paste0(wd, "TME_TIGIT/Seurat_subsets/Post_Annotation/scrna_immune_lymphoid_mb.rds"))

