library(Seurat)
library(readr)

Seurat_files = list.files("Seurat_objects")
Seurat_files[c(5,7,11)] = NA
Seurat_files = na.omit(Seurat_files)
Seurat_list = list()
for (i in Seurat_files) {
  load(paste0("Seurat_objects/", i))
  Seurat_list[[i]] = sce.m
  DefaultAssay(Seurat_list[[i]]) = "RNA"
  Seurat_list[[i]] = subset(Seurat_list[[i]], subset = percent_ribo<45) #After first run of analysis, lots of RP genes popped up, re-run with stricter RP filtering (10x genomics suggests that "normal" levels of ribosomial protein in T cells are 40-45%, using this as max value)
  Seurat_list[[i]] = SCTransform(Seurat_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
rm(sce.m)
features_sct = SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 5000)  
Seurat_list= PrepSCTIntegration(object.list = Seurat_list, anchor.features = features_sct)
anchors_sct <- FindIntegrationAnchors(object.list = Seurat_list, normalization.method = "SCT", 
                                             anchor.features = features_sct, dims = 1:30, reduction = "rpca", k.anchor = 20 )
rm(Seurat_list)
scrna <- IntegrateData(anchorset = anchors_sct, normalization.method = "SCT", dims = 1:30)
rm(anchors_sct)
scrna <- RunPCA(scrna, verbose = T, dims=1:30)
scrna = RunUMAP(scrna, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna = FindNeighbors(scrna, dims = 1:30)
i = seq(0,1,by = 0.1)
scrna= FindClusters(scrna, resolution = i, reduction = "umap")
write_rds(scrna, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/scrna_rpca.rds")