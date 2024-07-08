library(Seurat)
library(readr)
library(harmony)


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
integ_features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 3000) 
seurat_objects_first = Seurat_list[[1]]
Seurat_list[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(Seurat_list), merge.data = TRUE, project = "TME") 
rm(Seurat_list)
DefaultAssay(seurat_objects) = "SCT"
VariableFeatures(seurat_objects) = integ_features
scrna_merge = RunPCA (seurat_objects, verbose = FALSE, assay = "SCT", npcs= 50)
rm(seurat_objects)
scrna_merge = RunUMAP(scrna_merge, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna_merge = RunHarmony(scrna_merge, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_merge = FindNeighbors(object = scrna_merge, reduction = "harmony", dims = 1:30)
scrna_merge = RunUMAP(scrna_merge, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
i = seq(0.1, 1, by = 0.1)
scrna_merge = FindClusters(scrna_merge, resolution = i)
write_rds(scrna_merge, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/scrna_harmony.rds")