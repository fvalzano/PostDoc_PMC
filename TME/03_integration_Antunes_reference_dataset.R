library(Seurat)
library(readr)
library(dplyr)
library(SingleCellExperiment)
library(scater)


list = list.files("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM")
list[6] = NA
list = na.omit(list)
seurat_list = list()
seurat_objects = list()
Meta = list()
for (i in list) {
   seurat_list[[i]] = Read10X(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM/", i, "/filtered_feature_bc_matrix/filtered_feature_bc_matrices"))
   Meta[[i]] = read.csv(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM/", i, "/Meta.csv"))
   rownames(Meta[[i]]) = Meta[[i]]$cell
   seurat_objects[[i]] = CreateSeuratObject(counts = seurat_list[[i]], project = i, min.cells = 3, min.features = 200)
   seurat_objects[[i]] = seurat_objects[[i]][,colnames(seurat_objects[[i]]) %in% Meta[[i]]$cell]
}

for(i in names(seurat_objects)){
   ##### Normalize data
   seurat_objects[[i]] <- NormalizeData(seurat_objects[[i]],verbose = F)
   ##### HVG detection 
   seurat_objects[[i]] <- FindVariableFeatures(seurat_objects[[i]],verbose=T)
   ##### Scale data per gene
   seurat_objects[[i]] <- ScaleData(seurat_objects[[i]],verbose=T)
   ##### PCA
   seurat_objects[[i]] <- RunPCA(seurat_objects[[i]], features =VariableFeatures(seurat_objects[[i]]))
   seurat_objects[[i]] <- RunUMAP(seurat_objects[[i]], dims = 1:30, verbose=F)
   
   seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], metadata = Meta[[i]])
   seurat_objects[[i]]@meta.data = na.omit(seurat_objects[[i]]@meta.data)
}

for(i in names(seurat_objects)){
  write_rds(seurat_objects[[i]], paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Annotation_reference/Antunes_GBM/RDS_file/", i, ".rds"))
}

