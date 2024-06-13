library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(enrichplot)
library(enrichR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(scran)
#Load reference datasets and run SingleR autoannotation
reference = list.files("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM/RDS_file")
seurat_objects=list()

for (i in reference) {
  seurat_objects[[i]] = read_rds(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM/RDS_file/", i))
  Idents(seurat_objects[[i]]) = "cluster"
}
Ref1_Myeloid = subset(seurat_objects[[1]], idents = c("TAM 2", "Monocytes", "TAM 1", "prol. TAM", "DC"))
Ref2_Myeloid = subset(seurat_objects[[2]], idents = c("DC 1", "DC 2", "DC 3", "DC 4", "Monocytes", "prol. TAM", "TAM 1", "TAM 2"))
Ref3_Myeloid = subset(seurat_objects[[3]], idents = c("Hypoxic Mg-TAM", "IFN Mg-TAM", "Mg-TAM", "Phago/Lipid Mg-TAM"))
Ref1_Lymphoid = subset(seurat_objects[[1]], idents = c("B cells", "NK cells", "T cells"))
Ref2_Lymphoid = subset(seurat_objects[[2]], idents = c("B cells", "NK cells", "Plasma B", "Regulatory T cells", "T cells"))

scrna_immune_myeloid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune_myeloid_SCT.rds")
annotations_top = SingleR(test = scrna_immune_myeloid@assays$RNA$counts,
                          ref = list(Ref1_Myeloid@assays$RNA$counts,
                                     Ref2_Myeloid@assays$RNA$counts,
                                     Ref3_Myeloid@assays$RNA$counts),
                          labels = list(Ref1_Myeloid$cluster,
                                        Ref2_Myeloid$cluster,
                                        Ref3_Myeloid$cluster),
                          de.method="wilcox")
transfer.anno = as.data.frame(annotations_top$labels, row.names = rownames(annotations_top))
transfer.anno$`annotations_top$labels` = as.factor(transfer.anno$`annotations_top$labels`)
scrna_immune_myeloid <- AddMetaData(scrna_immune_myeloid, transfer.anno, col.name = "Antunes_annotations_top")
write_rds(scrna_immune_myeloid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_SingleR/scrna_immune_myeloid.rds")

scrna_immune_lymphoid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune_lymphoid_SCT.rds")
annotations_top = SingleR(test = scrna_immune_lymphoid@assays$RNA$counts,
                          ref = list(Ref1_Lymphoid@assays$RNA$counts,
                                     Ref2_Lymphoid@assays$RNA$counts),
                          labels = list(Ref1_Lymphoid$cluster,
                                        Ref2_Lymphoid$cluster),
                          de.method="wilcox")
transfer.anno = as.data.frame(annotations_top$labels, row.names = rownames(annotations_top))
transfer.anno$`annotations_top$labels` = as.factor(transfer.anno$`annotations_top$labels`)
scrna_immune_lymphoid <- AddMetaData(scrna_immune_lymphoid, transfer.anno, col.name = "Antunes_annotations_top")
write_rds(scrna_immune_lymphoid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_SingleR/scrna_immune_lymphoid.rds")

#Annotation_myeloid_v1
scrna_immune_myeloid = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/scrna_immune_myeloid.rds")
Idents(scrna_immune_myeloid) = "SCT_snn_res.0.6"
scrna_immune_myeloid = RenameIdents(scrna_immune_myeloid, c("0" = "TAM1",
                                                                     "1" = "TAM2",
                                                                     "2" = "TAM3",
                                                                     "3" = "TAM4",
                                                                     "4" = "Monocytes",
                                                                     "5" = "5",
                                                                     "6" = "DC",
                                                                     "7" = "TAM5",
                                                                     "8" = "TAM6",
                                                                     "9" = "TAM7",
                                                                     "10" = "TAM8",
                                                                     "11" = "prolif. TAM",
                                                                     "12" = "TAM9",
                                                                     "13" = "TAM10",
                                                                     "14" = "14",
                                                                     "15" = "15",
                                                                     "16" = "16",
                                                                     "17" = "17"),
                                                                     "18" = "18")
scrna_immune_myeloid$annotation_fv_v1 = scrna_immune_myeloid@active.ident
scrna_immune_myeloid = PrepSCTFindMarkers(scrna_immune_myeloid, assay="SCT")
DEG_annotation_fv_v1 = FindAllMarkers(scrna_immune_myeloid, 
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.15)
DEG_annotation_fv_v1 = DEG_annotation_fv_v1[DEG_annotation_fv_v1$p_val_adj <= 0.05,]
write.csv2(DEG_annotation_fv_v1, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/DEG/DEG_annotation_myeloid_fv_v1/Annotation_fv_v1.csv")
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()

#Split master DEG analysis per cluster and run enrichment analysis
for(i in unique(DEG_annotation_fv_v1$cluster)){
  DEG_list[[i]] = DEG_annotation_fv_v1[DEG_annotation_fv_v1$cluster ==i,]
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
  write.csv2(gse[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/EnrichmentResults/Enrichment_Myeloid_clusters/annotation_fv_v1/", i, ".csv"))
}
#Annotation_myeloid_v2
Idents(scrna_immune_myeloid) = "annotation_fv_v1"
scrna_immune_myeloid = RenameIdents(scrna_immune_myeloid, c("TAM1" = "Proinflammatory Microglia1", #CCL3,CCL4
                                                                     "TAM2" = "Immunosuppressive TAM1", #C1Q complex
                                                                     "TAM3" = "Immunosuppressive TAM2", #C1Q complex, TREM2
                                                                     "TAM4" = "Hypoxic TAM", #HIF1A
                                                                     "Monocytes" = "Monocytes",
                                                                     "TAM5" = "Immunosuppressive TAM3", #Not really displaying a strong fingerprint, by the disposition of the cells looks similar to immunosuppressive TAM
                                                                     "5" = "Astrocytes1",
                                                                     "DC" = "DC",
                                                                     "TAM6" = "Proinflammatory Microglia2", #CCL3,CCL4
                                                                     "TAM7" = "Immunosuppressive TAM3", #C1Q complex
                                                                     "prolif. TAM" = "prolif. TAM",
                                                                     "TAM8" = "HSP-TAM", #HSP genes, maybe remove?
                                                                     "TAM9" = "Activated TAM", #Dig deeper in here
                                                                     "TAM10" = "IFN-activated TAM",
                                                                     "15" = "Astrocytes2", 
                                                                     "14" = "Astrocytes3",
                                                                     "15" = "Astrocytes4",
                                                                     "16" = "Astrocytes5",
                                                                     "17" = "Astrocytes6",
                                                                     "18" = "Damaged Cells")) #MT genes scoring highest
scrna_immune_myeloid$annotation_fv_v2 = scrna_immune_myeloid@active.ident
DEG_annotation_fv_v2 = FindAllMarkers(scrna_immune_myeloid, 
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.15)
DEG_annotation_fv_v2 = DEG_annotation_fv_v2[DEG_annotation_fv_v2$p_val_adj <= 0.05,]
write.csv2(DEG_annotation_fv_v2, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/DEG/DEG_annotation_myeloid_fv_v2/Annotation_fv_v2.csv")

#Split master DEG analysis per cluster and run enrichment analysis
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
for(i in unique(DEG_annotation_fv_v2$cluster)){
  DEG_list[[i]] = DEG_annotation_fv_v2[DEG_annotation_fv_v2$cluster ==i,]
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
  write.csv2(gse[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/EnrichmentResults/Enrichment_Myeloid_clusters/annotation_fv_v2/", i, ".csv"))
}

write_rds(scrna_immune_myeloid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Annotation/scrna_immune_myeloid.rds")

#Annotation_lymphoid_v1
scrna_immune_lymphoid = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/scrna_immune_lymphoid.rds")
Idents(scrna_immune_lymphoid) = "SCT_snn_res.1.4"
scrna_immune_lymphoid = RenameIdents(scrna_immune_lymphoid, c("0" = "T Cells1",
                                                              "1" = "T Cells2",
                                                              "2" = "T Cells3",
                                                              "3" = "T Cells4",
                                                              "4" = "NK Cells",
                                                              "5" = "Intermediate Cells",
                                                              "6" = "T Cells6",
                                                              "7" = "7",
                                                              "8" = "T Cells7",
                                                              "9" = "B Cells",
                                                              "10" = "prolif. T Cells",
                                                              "11" = "Treg Cells",
                                                              "12" = "Plasma B Cells"))
scrna_immune_lymphoid$annotation_fv_v1 = scrna_immune_lymphoid@active.ident
scrna_immune_lymphoid = PrepSCTFindMarkers(scrna_immune_lymphoid, assay="SCT")
DEG_annotation_fv_v1 = FindAllMarkers(scrna_immune_lymphoid, 
                                      only.pos = T, 
                                      min.pct = 0.1,
                                      min.diff.pct = 0.1)
DEG_annotation_fv_v1 = DEG_annotation_fv_v1[DEG_annotation_fv_v1$p_val_adj <= 0.05,]
write.csv2(DEG_annotation_fv_v1, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/DEG/DEG_annotation_lymphoid_fv_v1/Annotation_fv_v1.csv")
#Split master DEG analysis per cluster and run enrichment analysis
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
for(i in unique(DEG_annotation_fv_v1$cluster)){
  DEG_list[[i]] = DEG_annotation_fv_v1[DEG_annotation_fv_v1$cluster ==i,]
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
  write.csv2(gse[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/EnrichmentResults/Enrichment_Lymphoid_clusters/annotation_fv_v1/", i, ".csv"))
}

#Annotation_lymphoid_v2
Idents(scrna_immune_lymphoid) = "annotation_fv_v1"
scrna_immune_lymphoid = RenameIdents(scrna_immune_lymphoid, c("T Cells1" = "Cytotoxic T Cells",  #GZM genes, CD8A
                                                              "T Cells2" = "Cytotoxic T Cells",
                                                              "T Cells3" = "Th Cells",   #JUN, JUNB -> Probably CD4 Th1 cells
                                                              "T Cells4" = "prolif. T Cells",
                                                              "NK Cells" = "NK Cells",
                                                              "Intermediate Cells" = "Intermediate Cells",  #Probably another name
                                                              "T Cells6" = "Naive T cells",   #Naive T cells according to AF
                                                              "7" = "RP-T Cells",   #Perhaps damaged cells? or Metabolic active?
                                                              "T Cells7" = "Dysfunctional T Cells",  #AF called them Artifacts - Are they though? Lots of genes involved in alternative RNA splicing
                                                              "B Cells" = "B Cells",
                                                              "prolif. T Cells" = "prolif. T Cells",
                                                              "NK Cells" = "NK Cells",
                                                              "Treg Cells" = "Treg Cells",
                                                              "Plasma B Cells" = "Plasma B Cells")) 
scrna_immune_lymphoid$annotation_fv_v2 = scrna_immune_lymphoid@active.ident
DEG_annotation_fv_v2 = FindAllMarkers(scrna_immune_lymphoid,
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.10)
DEG_annotation_fv_v2 = DEG_annotation_fv_v2[DEG_annotation_fv_v2$p_val_adj <= 0.05,]
write.csv2(DEG_annotation_fv_v2, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/DEG/DEG_annotation_lymphoid_fv_v2/Annotation_fv_v2.csv")
#Split master DEG analysis per cluster and run enrichment analysis
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
for(i in unique(DEG_annotation_fv_v2$cluster)){
  DEG_list[[i]] = DEG_annotation_fv_v2[DEG_annotation_fv_v2$cluster ==i,]
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
  write.csv2(gse[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/EnrichmentResults/Enrichment_Lymphoid_clusters/annotation_fv_v2/", i, ".csv"))
}
write_rds(scrna_immune_lymphoid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Annotation/scrna_immune_lymphoid.rds")

scrna = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/scrna_harmony.rds")
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/FullDatasetPlots.pdf", height = 5, width = 7.5)
DimPlot(scrna, reduction = "umap_harmony", label = T, group.by = "SCT_snn_res.0.4")
FeaturePlot(scrna, reduction = "umap_harmony", features="PTPRC", order = T, cols = c("#eeeeee", "#750000"))
dev.off()
scrna_immune_myeloid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_SingleR/scrna_immune_myeloid.rds")
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/MyeloidPlots.pdf", height = 5, width = 7.5)
DimPlot(scrna_immune_myeloid, reduction = "umap_harmony", label = T, group.by = "SCT_snn_res.0.6")
DimPlot(scrna_immune_myeloid, reduction = "umap_harmony", label = F, group.by = "Antunes_annotations_top")
DimPlot(scrna_immune_myeloid, reduction = "umap_harmony", label = F, group.by = "annotation_fv_v2")
FeaturePlot(scrna_immune_myeloid, reduction = "umap_harmony", features="CCL4", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_myeloid, reduction = "umap_harmony", features="C1QA", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_myeloid, reduction = "umap_harmony", features="TREM2", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_myeloid, reduction = "umap_harmony", features="GFAP", order = T, cols = c("#eeeeee", "#750000"))
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/MyeloidPlotsSplit.pdf", height = 5, width = 12.5)
DimPlot(scrna_immune_myeloid, reduction = "umap_harmony", label = F, group.by = "annotation_fv_v2", split.by = "Entity")
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/EnrichmentPlots_Myeloid.pdf", height = 5, width = 5)
for (i in names(gse)) {
  p = dotplot(gse[[i]], showCategory=10)
  plot(p)
}
dev.off()
scrna_immune_lymphoid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_SingleR/scrna_immune_lymphoid.rds")
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/LymphoidPlots.pdf", height = 5, width = 7.5)
DimPlot(scrna_immune_lymphoid, reduction = "umap_harmony", label = T, group.by = "SCT_snn_res.1.4")
DimPlot(scrna_immune_lymphoid, reduction = "umap_harmony", label = F, group.by = "Antunes_annotations_top")
DimPlot(scrna_immune_lymphoid, reduction = "umap_harmony", label = F, group.by = "annotation_fv_v2")
FeaturePlot(scrna_immune_lymphoid, reduction = "umap_harmony", features="CD3D", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_lymphoid, reduction = "umap_harmony", features="CD8A", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_lymphoid, reduction = "umap_harmony", features="GZMK", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_lymphoid, reduction = "umap_harmony", features="CD79A", order = T, cols = c("#eeeeee", "#750000"))
FeaturePlot(scrna_immune_lymphoid, reduction = "umap_harmony", features="NKG7", order = T, cols = c("#eeeeee", "#750000"))
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/LymphoidPlotsSplit.pdf", height = 5, width = 12.5)
DimPlot(scrna_immune_lymphoid, reduction = "umap_harmony", label = F, group.by = "annotation_fv_v2", split.by = "Entity")
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Plots/EnrichmentPlots_Lymphoid.pdf", height = 5, width = 5)
for (i in names(gse)) {
  p = dotplot(gse[[i]], showCategory=10)
  plot(p)
}
dev.off()
