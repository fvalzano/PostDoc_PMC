library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)

#Subcluster full dataset in Immune subclusters
scrna = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/scrna_harmony.rds")
Idents(scrna) = "SCT_snn_res.0.4"
scrna_immune = subset(scrna, idents = c("0", "28", "41"))
rm(scrna)
scrna_immune_list = SplitObject(scrna_immune, split.by = "Dataset")
hvg_immune = SelectIntegrationFeatures(scrna_immune_list, nfeatures = 3000)
rm(scrna_immune_list)
scrna_immune = RunPCA (scrna_immune, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_immune)
scrna_immune = RunUMAP(scrna_immune, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_immune = FindNeighbors(object = scrna_immune, reduction = "harmony", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_immune = FindClusters(scrna_immune, resolution = i)
write_rds(scrna_immune, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune.rds")

#Myeloid compartment
scrna_immune = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/scrna_immune.rds")
Idents(scrna_immune) = "SCT_snn_res.0.2"
scrna_immune_myeloid = subset(scrna_immune, idents = c("0", "2", "3", "4", "5", "6", "7"))
scrna_immune_myeloid_list = SplitObject(scrna_immune_myeloid, split.by = "Dataset")
for (i in names(scrna_immune_myeloid_list)) {
  DefaultAssay(scrna_immune_myeloid_list[[i]]) = "RNA"
  scrna_immune_myeloid_list[[i]] = SCTransform(scrna_immune_myeloid_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_immune_myeloid = SelectIntegrationFeatures(scrna_immune_myeloid_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_immune_myeloid, value = TRUE)
subtract<-which(hvg_immune_myeloid %in% ribo.genes)
hvg_immune_myeloid_filtered<-hvg_immune_myeloid[-subtract]
scrna_immune_myeloid = merge(x = scrna_immune_myeloid_list[[1]], y= scrna_immune_myeloid_list[-1], merge.data = TRUE, project = "TME_myeloid") 
#Harmony
scrna_immune_myeloid = RunPCA (scrna_immune_myeloid, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_immune_myeloid_filtered)
scrna_immune_myeloid = RunHarmony(scrna_immune_myeloid, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_immune_myeloid = RunUMAP(scrna_immune_myeloid, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_immune_myeloid = FindNeighbors(object = scrna_immune_myeloid, reduction = "harmony", dims = 1:30)
i = seq(0.1, 1, by = 0.1)
scrna_immune_myeloid = FindClusters(scrna_immune_myeloid, resolution = i)
#NMF
#scrna_immune_myeloid = runNMF(scrna_immune_myeloid, k = 30, hvg = hvg_immune_myeloid_filtered, assay = "SCT", reduction="nmf")
#scrna_immune_myeloid = RunUMAP(scrna_immune_myeloid, reduction = "NMF", dims=1:30, reduction.name = "umap_nmf", reduction.key = "nmfUMAP_")
#scrna_immune_myeloid = FindNeighbors(object = scrna_immune_myeloid, reduction = "nmf", dims = 1:30)
#scrna_immune_myeloid_list = SplitObject(scrna_immune_myeloid, split.by = "Dataset")
#hvg_immune_myeloid = SelectIntegrationFeatures(scrna_immune_myeloid_list, nfeatures = 3000)
#ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_immune_ymphoid, value = TRUE)
#subtract<-which(hvg_immune_ymphoid %in% ribo.genes)
#hvg_immune_myeloid_filtered<-hvg_immune_myeloid[-subtract]
#scrna_immune_myeloid = runNMF(scrna_immune_myeloid, k = 30, hvg = hvg_immune_myeloid_filtered, assay = "SCT")
#scrna_immune_myeloid@reductions$NMF
#scrna_immune_myeloid = RunUMAP(scrna_immune_myeloid, reduction = "NMF", dims=1:30, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")
#geneNMF.programs <- multiNMF(scrna_immune_myeloid_list, assay="SCT", slot="data", k=1:30, L1=c(0,0),
#                    do_centering=TRUE, nfeatures = 3000)
#geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                        nprograms=15,
#                                        max.genes=50,
#                                        hclust.method="ward.D2",
#                                        min.confidence=0.3)
#mp.genes <- geneNMF.metaprograms$metaprograms.genes
#scrna_immune_myeloid <- AddModuleScore_UCell(scrna_immune_myeloid, features = mp.genes, assay="SCT", ncores=4, name = "")
#matrix <- scrna_immune_myeloid@meta.data[,names(mp.genes)]
#dimred <- as.matrix(matrix)
#colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
#scrna_immune_myeloid@reductions[["NMF_MPs"]] <- new("DimReduc",
#                                         cell.embeddings = dimred,
#                                         assay.used = "RNA",
#                                         key = "MP_",
#                                         global = FALSE)
#scrna_immune_myeloid <- RunUMAP(scrna_immune_myeloid, reduction="NMF_MPs", dims=1:length(scrna_immune_myeloid@reductions[["NMF_MPs"]]),
#              metric = "euclidean", reduction.name = "umap_MP")

write_rds(scrna_immune_myeloid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune_myeloid_SCT.rds")

scrna_immune = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune.rds")
Idents(scrna_immune) = "SCT_snn_res.0.2"
scrna_immune_myeloid = subset(scrna_immune, idents = c("0", "2", "3", "4", "5", "6", "7"))
DefaultAssay(scrna_immune_myeloid) = "RNA"
FindVariableFeatures(scrna_immune_myeloid)
scrna_immune_myeloid = NormalizeData(scrna_immune_myeloid, normalization.method="LogNormalize", scale.factor=10000)
scrna_immune_myeloid = ScaleData(scrna_immune_myeloid, features=rownames(scrna_immune_myeloid), vars.to.regress = c("percent_mito", "percent_ribo"))
scrna_immune_myeloid = RunPCA (scrna_immune_myeloid, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_immune_myeloid = RunUMAP(scrna_immune_myeloid, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna_immune_myeloid = FindNeighbors(object = scrna_immune_myeloid, reduction = "pca", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
#scrna_immune_myeloid = FindClusters(scrna_immune_myeloid, resolution = i)
write_rds(scrna_immune_myeloid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_myeloid_standardwfl.rds")

#Lymphoid compartment
scrna_immune = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune.rds")
Idents(scrna_immune) = "SCT_snn_res.0.2"
scrna_immune_lymphoid = subset(scrna_immune, idents = c("1", "8"))
scrna_immune_lymphoid_list = SplitObject(scrna_immune_lymphoid, split.by = "Dataset")
for (i in names(scrna_immune_lymphoid_list)) {
  DefaultAssay(scrna_immune_lymphoid_list[[i]]) = "RNA"
  scrna_immune_lymphoid_list[[i]] = SCTransform(scrna_immune_lymphoid_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_immune_lymphoid = SelectIntegrationFeatures(scrna_immune_lymphoid_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_immune_lymphoid, value = TRUE)
subtract<-which(hvg_immune_lymphoid %in% ribo.genes)
hvg_immune_lymphoid_filtered<-hvg_immune_lymphoid[-subtract]
scrna_immune_lymphoid = merge(x = scrna_immune_lymphoid_list[[1]], y= scrna_immune_lymphoid_list[-1], merge.data = TRUE, project = "TME_lymphoid") 
#Harmony
scrna_immune_lymphoid = RunPCA (scrna_immune_lymphoid, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_immune_lymphoid_filtered)
scrna_immune_lymphoid = RunHarmony(scrna_immune_lymphoid, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_immune_lymphoid = RunUMAP(scrna_immune_lymphoid, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_immune_lymphoid = FindNeighbors(object = scrna_immune_lymphoid, reduction = "harmony", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
scrna_immune_lymphoid = FindClusters(scrna_immune_lymphoid, resolution = i)
#NMF
#scrna_immune_lymphoid = runNMF(scrna_immune_lymphoid, k = 30, hvg = hvg_immune_lymphoid_filtered, assay = "SCT", reduction="nmf")
#scrna_immune_lymphoid = RunUMAP(scrna_immune_lymphoid, reduction = "NMF", dims=1:30, reduction.name = "umap_nmf", reduction.key = "nmfUMAP_")
#scrna_immune_lymphoid = FindNeighbors(object = scrna_immune_lymphoid, reduction = "nmf", dims = 1:30)
#scrna_immune_lymphoid_list = SplitObject(scrna_immune_lymphoid, split.by = "Dataset")
#hvg_immune_lymphoid = SelectIntegrationFeatures(scrna_immune_lymphoid_list, nfeatures = 3000)
#ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_immune_ymphoid, value = TRUE)
#subtract<-which(hvg_immune_ymphoid %in% ribo.genes)
#hvg_immune_lymphoid_filtered<-hvg_immune_lymphoid[-subtract]
#scrna_immune_lymphoid = runNMF(scrna_immune_lymphoid, k = 30, hvg = hvg_immune_lymphoid_filtered, assay = "SCT")
#scrna_immune_lymphoid@reductions$NMF
#scrna_immune_lymphoid = RunUMAP(scrna_immune_lymphoid, reduction = "NMF", dims=1:30, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")
#geneNMF.programs <- multiNMF(scrna_immune_lymphoid_list, assay="SCT", slot="data", k=1:30, L1=c(0,0),
#                    do_centering=TRUE, nfeatures = 3000)
#geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                        nprograms=15,
#                                        max.genes=50,
#                                        hclust.method="ward.D2",
#                                        min.confidence=0.3)
#mp.genes <- geneNMF.metaprograms$metaprograms.genes
#scrna_immune_lymphoid <- AddModuleScore_UCell(scrna_immune_lymphoid, features = mp.genes, assay="SCT", ncores=4, name = "")
#matrix <- scrna_immune_lymphoid@meta.data[,names(mp.genes)]
#dimred <- as.matrix(matrix)
#colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
#scrna_immune_lymphoid@reductions[["NMF_MPs"]] <- new("DimReduc",
#                                         cell.embeddings = dimred,
#                                         assay.used = "RNA",
#                                         key = "MP_",
#                                         global = FALSE)
#scrna_immune_lymphoid <- RunUMAP(scrna_immune_lymphoid, reduction="NMF_MPs", dims=1:length(scrna_immune_lymphoid@reductions[["NMF_MPs"]]),
#              metric = "euclidean", reduction.name = "umap_MP")
write_rds(scrna_immune_lymphoid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune_lymphoid_SCT.rds")

scrna_immune = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune.rds")
Idents(scrna_immune) = "SCT_snn_res.0.2"
scrna_immune_lymphoid = subset(scrna_immune, idents = c("1", "8"))
DefaultAssay(scrna_immune_lymphoid) = "RNA"
FindVariableFeatures(scrna_immune_lymphoid)
scrna_immune_lymphoid = NormalizeData(scrna_immune_lymphoid, normalization.method="LogNormalize", scale.factor=10000)
scrna_immune_lymphoid = ScaleData(scrna_immune_lymphoid, features=rownames(scrna_immune_lymphoid), vars.to.regress = c("percent_mito", "percent_ribo"))
scrna_immune_lymphoid = RunPCA (scrna_immune_lymphoid, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_immune_lymphoid = RunUMAP(scrna_immune_lymphoid, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna_immune_lymphoid = FindNeighbors(object = scrna_immune_lymphoid, reduction = "pca", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
#scrna_immune_lymphoid = FindClusters(scrna_immune_lymphoid, resolution = i)
write_rds(scrna_immune_lymphoid, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_lymphoid_standardwfl.rds")
