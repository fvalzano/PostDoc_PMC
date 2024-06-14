library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)

###All
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_ep = read_rds(paste0(wd,"Seurat_subsets/Post_Entity_Splitting/scrna_harmony_Ependymoma.rds"))
scrna_ep_list = SplitObject(scrna_ep, split.by = "Dataset")
for (i in names(scrna_ep_list)) {
  DefaultAssay(scrna_ep_list[[i]]) = "RNA"
  scrna_ep_list[[i]] = SCTransform(scrna_ep_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg= SelectIntegrationFeatures(scrna_ep_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg, value = TRUE)
subtract<-which(hvg %in% ribo.genes)
hvg_filtered<-hvg[-subtract]
scrna_ep = merge(x = scrna_ep_list[[1]], y= scrna_ep_list[-1], merge.data = TRUE, project = "Ep") 
scrna_ep = RunPCA (scrna_ep, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_filtered)
scrna_ep = RunHarmony(scrna_ep, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_ep = RunUMAP(scrna_ep, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_ep = FindNeighbors(object = scrna_ep, reduction = "harmony", dims = 1:30)
i = seq(0.1, 1, by = 0.1)
scrna_ep = FindClusters(scrna_ep, resolution = i)
write_rds(scrna_ep, paste0(wd, "TME_B7H3/Seurat_subsets/scrna_ep.rds"))

###Myeloid
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_myeloid = read_rds(paste0(wd,"Seurat_subsets/Post_Annotation/scrna_immune_myeloid.rds"))
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
i = seq(0.1, 1, by = 0.1)
scrna_myeloid = FindClusters(scrna_myeloid, resolution = i)
write_rds(scrna_myeloid, paste0(wd, "TME_B7H3/Seurat_subsets/scrna_immune_myeloid.rds"))

###Lymphoid
wd = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/"
scrna_lymphoid = read_rds(paste0(wd,"Seurat_subsets/Post_Annotation/scrna_immune_lymphoid.rds"))
scrna_lymphoid_list = SplitObject(scrna_lymphoid, split.by = "Dataset")
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
scrna_lymphoid = FindNeighbors(object = scrna_lymphoid, reduction = "harmony", dims = 1:30)
i = seq(0.1, 1, by = 0.1)
scrna_lymphoid = FindClusters(scrna_lymphoid, resolution = i)
write_rds(scrna_lymphoid, paste0(wd, "TME_B7H3/Seurat_subsets/scrna_immune_lymphoid.rds"))

