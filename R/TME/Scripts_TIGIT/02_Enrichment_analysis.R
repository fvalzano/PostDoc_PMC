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

###Myeloid
wd="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/"
scrna_myeloid = readRDS(paste0(wd, "Seurat_subsets/Post_Annotation/scrna_immune_myeloid_mb.rds"))
Idents(scrna_myeloid) = "annotation_fv_v2"
scrna_myeloid = PrepSCTFindMarkers(scrna_myeloid)
DEG_annotation = FindAllMarkers(scrna_myeloid, 
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.15)
DEG_annotation = DEG_annotation[DEG_annotation$p_val_adj <= 0.05,]
write.csv2 = paste0(wd, "DEG/Myeloid/DEG_Clusters.csv")
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
for(i in unique(DEG_annotation$cluster)){
  DEG_list[[i]] = DEG_annotation[DEG_annotation$cluster ==i,]
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
  write.csv2(gse[[i]], paste0((wd), "Enrichment/Myeloid/Enrichment_Myeloid_clusters_", i, ".csv"))
}

###Lymphoid
wd="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/TME_TIGIT/"
scrna_lymphoid = readRDS(paste0(wd, "Seurat_subsets/Post_Annotation/scrna_immune_lymphoid_mb.rds"))
Idents(scrna_lymphoid) = "annotation_fv_v2"
scrna_lymphoid = PrepSCTFindMarkers(scrna_lymphoid)
DEG_annotation = FindAllMarkers(scrna_lymphoid, 
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.1)
DEG_annotation = DEG_annotation[DEG_annotation$p_val_adj <= 0.05,]
write.csv2 = paste0(wd, "DEG/Lymphoid/DEG_Clusters.csv")
DEG_list = list()
gse= list()
DEG_list_IDs = list()
DEG_vector = list()
for(i in unique(DEG_annotation$cluster)){
  DEG_list[[i]] = DEG_annotation[DEG_annotation$cluster ==i,]
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
  write.csv2(gse[[i]], paste0((wd), "Enrichment/Lymphoid/Enrichment_lymphoid_clusters_", i, ".csv"))
}
