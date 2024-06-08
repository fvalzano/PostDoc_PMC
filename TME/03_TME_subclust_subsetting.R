library(Seurat)
library(readr)
library(SingleR)
library(harmony)
library(GeneNMF)
library(UCell)

#Subcluster full dataset in subclusters of interest
scrna = read_rds("YourDirectory/scrna_object.rds")
Idents(scrna) = "SCT_snn_res.0.4"
scrna_subset = subset(scrna, idents = c("0", "28", "41"))
rm(scrna)
scrna_subset_list = SplitObject(scrna_subset, split.by = "Dataset")
hvg = SelectIntegrationFeatures(scrna_subset_list, nfeatures = 3000)
rm(scrna_subset_list)
scrna_subset = RunPCA (scrna_subset, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg)
scrna_subset = RunUMAP(scrna_subset, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_subset = FindNeighbors(object = scrna_subset, reduction = "harmony", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_subset = FindClusters(scrna_subset, resolution = i)
write_rds(scrna_subset, "YourDirectory/scrna_subset.rds")

#my compartment
scrna_subset = read_rds("YourDirectory/scrna_subset.rds")
Idents(scrna_subset) = "SCT_snn_res.0.2"
scrna_subset_my = subset(scrna_subset, idents = c("0", "2", "3", "4", "5", "6", "7"))
scrna_subset_my_list = SplitObject(scrna_subset_my, split.by = "Dataset")
for (i in names(scrna_subset_my_list)) {
  DefaultAssay(scrna_subset_my_list[[i]]) = "RNA"
  scrna_subset_my_list[[i]] = SCTransform(scrna_subset_my_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_my = SelectIntegrationFeatures(scrna_subset_my_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_my, value = TRUE)
subtract<-which(hvg_my %in% ribo.genes)
hvg_my_filtered<-hvg_my[-subtract]
scrna_subset_my = merge(x = scrna_subset_my_list[[1]], y= scrna_subset_my_list[-1], merge.data = TRUE, project = "TME_my") 
#Harmony
scrna_subset_my = RunPCA (scrna_subset_my, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_my_filtered)
scrna_subset_my = RunHarmony(scrna_subset_my, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_subset_my = RunUMAP(scrna_subset_my, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_subset_my = FindNeighbors(object = scrna_subset_my, reduction = "harmony", dims = 1:30)
i = seq(0.1, 1, by = 0.1)
scrna_subset_my = FindClusters(scrna_subset_my, resolution = i)
#NMF
#scrna_subset_my = runNMF(scrna_subset_my, k = 30, hvg = hvg_my_filtered, assay = "SCT", reduction="nmf")
#scrna_subset_my = RunUMAP(scrna_subset_my, reduction = "NMF", dims=1:30, reduction.name = "umap_nmf", reduction.key = "nmfUMAP_")
#scrna_subset_my = FindNeighbors(object = scrna_subset_my, reduction = "nmf", dims = 1:30)
#scrna_subset_my_list = SplitObject(scrna_subset_my, split.by = "Dataset")
#hvg_my = SelectIntegrationFeatures(scrna_subset_my_list, nfeatures = 3000)
#ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_ymphoid, value = TRUE)
#subtract<-which(hvg_ymphoid %in% ribo.genes)
#hvg_my_filtered<-hvg_my[-subtract]
#scrna_subset_my = runNMF(scrna_subset_my, k = 30, hvg = hvg_my_filtered, assay = "SCT")
#scrna_subset_my@reductions$NMF
#scrna_subset_my = RunUMAP(scrna_subset_my, reduction = "NMF", dims=1:30, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")
#geneNMF.programs <- multiNMF(scrna_subset_my_list, assay="SCT", slot="data", k=1:30, L1=c(0,0),
#                    do_centering=TRUE, nfeatures = 3000)
#geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                        nprograms=15,
#                                        max.genes=50,
#                                        hclust.method="ward.D2",
#                                        min.confidence=0.3)
#mp.genes <- geneNMF.metaprograms$metaprograms.genes
#scrna_subset_my <- AddModuleScore_UCell(scrna_subset_my, features = mp.genes, assay="SCT", ncores=4, name = "")
#matrix <- scrna_subset_my@meta.data[,names(mp.genes)]
#dimred <- as.matrix(matrix)
#colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
#scrna_subset_my@reductions[["NMF_MPs"]] <- new("DimReduc",
#                                         cell.embeddings = dimred,
#                                         assay.used = "RNA",
#                                         key = "MP_",
#                                         global = FALSE)
#scrna_subset_my <- RunUMAP(scrna_subset_my, reduction="NMF_MPs", dims=1:length(scrna_subset_my@reductions[["NMF_MPs"]]),
#              metric = "euclidean", reduction.name = "umap_MP")

write_rds(scrna_subset_my, "YourDirectory/scrna_subset_my_SCT.rds")

#StandardWorkflow
scrna_subset = read_rds("YourDirectory/scrna_subset.rds")
Idents(scrna_subset) = "SCT_snn_res.0.2"
scrna_subset_my = subset(scrna_subset, idents = c("0", "2", "3", "4", "5", "6", "7"))
DefaultAssay(scrna_subset_my) = "RNA"
FindVariableFeatures(scrna_subset_my)
scrna_subset_my = NormalizeData(scrna_subset_my, normalization.method="LogNormalize", scale.factor=10000)
scrna_subset_my = ScaleData(scrna_subset_my, features=rownames(scrna_subset_my), vars.to.regress = c("percent_mito", "percent_ribo"))
scrna_subset_my = RunPCA (scrna_subset_my, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_subset_my = RunUMAP(scrna_subset_my, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna_subset_my = FindNeighbors(object = scrna_subset_my, reduction = "pca", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
#scrna_subset_my = FindClusters(scrna_subset_my, resolution = i)
write_rds(scrna_subset_my, "YourDirectory/scrna_my_standardwfl.rds")

#ly compartment
scrna_subset = read_rds("YourDirectory/scrna_subset.rds")
Idents(scrna_subset) = "SCT_snn_res.0.2"
scrna_subset_ly = subset(scrna_subset, idents = c("1", "8"))
scrna_subset_ly_list = SplitObject(scrna_subset_ly, split.by = "Dataset")
for (i in names(scrna_subset_ly_list)) {
  DefaultAssay(scrna_subset_ly_list[[i]]) = "RNA"
  scrna_subset_ly_list[[i]] = SCTransform(scrna_subset_ly_list[[i]], vars.to.regress = c("percent_mito", "percent_ribo"), vst.flavor = "v2", assay = "RNA")
}
hvg_ly = SelectIntegrationFeatures(scrna_subset_ly_list, nfeatures = 3000)
ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_ly, value = TRUE)
subtract<-which(hvg_ly %in% ribo.genes)
hvg_ly_filtered<-hvg_ly[-subtract]
scrna_subset_ly = merge(x = scrna_subset_ly_list[[1]], y= scrna_subset_ly_list[-1], merge.data = TRUE, project = "TME_ly") 
#Harmony
scrna_subset_ly = RunPCA (scrna_subset_ly, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg_ly_filtered)
scrna_subset_ly = RunHarmony(scrna_subset_ly, group.by.vars = "Dataset", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_subset_ly = RunUMAP(scrna_subset_ly, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_subset_ly = FindNeighbors(object = scrna_subset_ly, reduction = "harmony", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
scrna_subset_ly = FindClusters(scrna_subset_ly, resolution = i)
#NMF
#scrna_subset_ly = runNMF(scrna_subset_ly, k = 30, hvg = hvg_ly_filtered, assay = "SCT", reduction="nmf")
#scrna_subset_ly = RunUMAP(scrna_subset_ly, reduction = "NMF", dims=1:30, reduction.name = "umap_nmf", reduction.key = "nmfUMAP_")
#scrna_subset_ly = FindNeighbors(object = scrna_subset_ly, reduction = "nmf", dims = 1:30)
#scrna_subset_ly_list = SplitObject(scrna_subset_ly, split.by = "Dataset")
#hvg_ly = SelectIntegrationFeatures(scrna_subset_ly_list, nfeatures = 3000)
#ribo.genes <- grep(pattern = "^RP[SL]", x = hvg_ymphoid, value = TRUE)
#subtract<-which(hvg_ymphoid %in% ribo.genes)
#hvg_ly_filtered<-hvg_ly[-subtract]
#scrna_subset_ly = runNMF(scrna_subset_ly, k = 30, hvg = hvg_ly_filtered, assay = "SCT")
#scrna_subset_ly@reductions$NMF
#scrna_subset_ly = RunUMAP(scrna_subset_ly, reduction = "NMF", dims=1:30, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")
#geneNMF.programs <- multiNMF(scrna_subset_ly_list, assay="SCT", slot="data", k=1:30, L1=c(0,0),
#                    do_centering=TRUE, nfeatures = 3000)
#geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
#                                        nprograms=15,
#                                        max.genes=50,
#                                        hclust.method="ward.D2",
#                                        min.confidence=0.3)
#mp.genes <- geneNMF.metaprograms$metaprograms.genes
#scrna_subset_ly <- AddModuleScore_UCell(scrna_subset_ly, features = mp.genes, assay="SCT", ncores=4, name = "")
#matrix <- scrna_subset_ly@meta.data[,names(mp.genes)]
#dimred <- as.matrix(matrix)
#colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
#scrna_subset_ly@reductions[["NMF_MPs"]] <- new("DimReduc",
#                                         cell.embeddings = dimred,
#                                         assay.used = "RNA",
#                                         key = "MP_",
#                                         global = FALSE)
#scrna_subset_ly <- RunUMAP(scrna_subset_ly, reduction="NMF_MPs", dims=1:length(scrna_subset_ly@reductions[["NMF_MPs"]]),
#              metric = "euclidean", reduction.name = "umap_MP")
write_rds(scrna_subset_ly, "YourDirectory/scrna_subset_ly_SCT.rds")

scrna_subset = read_rds("YourDirectory/scrna_subset.rds")
Idents(scrna_subset) = "SCT_snn_res.0.2"
scrna_subset_ly = subset(scrna_subset, idents = c("1", "8"))
DefaultAssay(scrna_subset_ly) = "RNA"
FindVariableFeatures(scrna_subset_ly)
scrna_subset_ly = NormalizeData(scrna_subset_ly, normalization.method="LogNormalize", scale.factor=10000)
scrna_subset_ly = ScaleData(scrna_subset_ly, features=rownames(scrna_subset_ly), vars.to.regress = c("percent_mito", "percent_ribo"))
scrna_subset_ly = RunPCA (scrna_subset_ly, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_subset_ly = RunUMAP(scrna_subset_ly, reduction = "pca", dims = 1:30, reduction.name = "umap")
scrna_subset_ly = FindNeighbors(object = scrna_subset_ly, reduction = "pca", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
#scrna_subset_ly = FindClusters(scrna_subset_ly, resolution = i)
write_rds("YourDirectory/scrna_ly_standardwfl.rds")
