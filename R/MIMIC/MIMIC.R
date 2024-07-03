---
title: "MIMIC_MB"
author: "Francesco Valzano"
date: "1/4/2024"
output: html_document
---
#All
##Preprocess
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r libraries, include=FALSE}
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(readr)
library(SingleR)
library(pheatmap)
library(dittoSeq)
library(SeuratWrappers)
library(monocle3)
library(DoubletFinder)
library(magrittr)
library(tradeSeq)
library(DropletUtils)
library(enrichR)
library(infercnv)
library(tidyr)
library(str2str)
library(NGCHM)
library(infercnvNGCHM)
library(futile.logger)
library(scater)
```
```{r MIMIC input files, include=FALSE}
base_directory = getwd() #Make sure you have the file you want to analyze in this folder: Name the folder with sequencing facility ID and put filtered_feature_bc_matrix in it, the script will fetch the data that it needs
filenames = c(paste0("LX", 430:441), paste0("LX",500:507))
raw_reads = list()
seurat_objects=list()


for (i in filenames) {
  # Create the file name for each sample
  file_dir = paste0(base_directory,"/Rstudio_Test1/MIMIC/Data/10x_runs/", i, "/filtered_feature_bc_matrix")
  # Read the sample data from the file using Read10x
  sample_data = Read10X(file_dir)
  raw_reads[[i]] = sample_data
  seurat_objects[[i]] = CreateSeuratObject(counts = raw_reads[[i]], project = i, min.cells = 3, min.features = 150)
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[rb.genes,])/Matrix::colSums(Assay)*100
  seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
  seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
  }
rm(sample_data, raw_reads) #Housekeeping

```
```{r QC, include=FALSE}
#Save QC plots
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/QC.pdf", width=5, height=5)
for (i in filenames) {
  #calculation of threshold is performed either by IsOutlier or via expliciting the formula
  min.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) - 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) + 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nCount.thr = median(seurat_objects[[i]]$nCount_RNA) + 3*mad(seurat_objects[[i]]$nCount_RNA)
  QC = seurat_objects[[i]]@meta.data
  QC$nCount_RNA = QC$nCount_RNA
  QC$nFeature_RNA = QC$nFeature_RNA
  QC$Dropouts = ifelse(QC$nCount_RNA<max.nCount.thr&
                         QC$nFeature_RNA>min.nFeature.thr&
                         QC$nFeature_RNA<max.nFeature.thr, FALSE, TRUE)
  p = QC %>%
    arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, colour=percent.mt, shape = Dropouts)) + 
    geom_point() + 
    scale_shape_manual(values = c(16,17), 
                       labels = c(paste0("Retained (",sum(QC$Dropouts == FALSE), " cells)"),
                                  paste0("Dropouts (",sum(QC$Dropouts == TRUE), " cells)"))) +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle(paste0("QC metrics ", unique(QC$orig.ident)))+
    geom_vline(xintercept = max.nCount.thr)+
    geom_hline(yintercept = c(min.nFeature.thr,max.nFeature.thr))+
    ylim(min.nFeature.thr-1,max(QC$nFeature_RNA))
  print(p)
}
dev.off()
rm(p,QC)
#Apply QC
qc_seurat_objects <- list()
for (i in filenames) {
  seurat_objects[[i]]$nFeature.thr = isOutlier(seurat_objects[[i]]$nFeature_RNA, nmads=3, type="both", log=FALSE)
   seurat_objects[[i]]$nCount.thr = isOutlier(seurat_objects[[i]]$nCount_RNA, nmads=3, type="higher", log=FALSE)
   seurat_objects[[i]] <- seurat_objects[[i]][,!(seurat_objects[[i]]$nCount.thr |  seurat_objects[[i]]$nFeature.thr)]
   seurat_objects[[i]] = subset(seurat_objects[[i]], subset = percent.mt <5 &
                                  percent.ribo<5)
  qc_seurat_objects[[i]] = subset(seurat_objects[[i]], subset = nFeature_RNA>min.nFeature.thr & 
                                 nFeature_RNA<max.nFeature.thr & 
                                 nCount_RNA<max.nCount.thr & 
                                 percent.mt<5 &
                                 percent.ribo<5)
}
rm(seurat_objects)


```
##Strategies
```{r Data Merging -Important in Seurat v5- SCT Normalization, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects =  merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
seurat_objects= CellCycleScoring(seurat_objects, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_mimic<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "pca", dims = 1:30)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "pca", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
```
```{r Harmony Integration}
scrna_mimic <- IntegrateLayers(object = scrna_mimic, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
DimPlot(scrna_mimic, group.by = "orig.ident")
write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
```
#Decided for SCT normalization and Harmony Integration, rest of the strategies are blanked from the main code
#```{r Data Merging -Important in Seurat v5- SCT Normalization + CCA Integration, eval=FALSE, include=FALSE}
#seurat_objects=qc_seurat_objects
#seurat_objects_first = seurat_objects[[1]]
#seurat_objects[[1]] = NULL
#seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
#seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2")
#seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
#scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "SCT")
#scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:30)
#i = seq(0.2, 1, by = 0.2)
#scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
#scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:30, reduction.name = "umap")
##write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_SCT_CCA/scrna.rds")
#```
#```{r Data Merging -Important in Seurat v5- Log Normalization + Harmony Integration, eval=FALSE, include=FALSE}
#seurat_objects=qc_seurat_objects
#seurat_objects_first = seurat_objects[[1]]
#seurat_objects[[1]] = NULL
#seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
#seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
#seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects))#, vars.to.regress = c("percent.mt", "percent.ribo"))
#seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
#scrna_mimic <- IntegrateLayers(object = seurat_objects, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
#scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
#i = seq(0.2, 1, by = 0.2)
#scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
#scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap")
##write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_Harmony/scrna.rds")
#DimPlot(scrna_mimic, split.by = "orig.ident")
#```
#```{r Data Merging -Important in Seurat v5- Log Normalization + CCA Integration, eval=FALSE, include=FALSE}
#seurat_objects=qc_seurat_objects
#seurat_objects_first = seurat_objects[[1]]
#seurat_objects[[1]] = NULL
#seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
#seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
#seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects), vars.to.regress = c("percent.mt", "percent.ribo"))
#seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
#scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
#scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:50)
#i = seq(0.2, 1, by = 0.2)
#scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
#scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:50, reduction.name = "umap")
##write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_CCA/scrna.rds")
#DimPlot(scrna_mimic, group.by = "orig.ident")
#```
##Metadata curation
```{r MimicID}
metadata_sequencing = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 1)
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata = list()
for (i in filenames) {
  metadata[[i]]= metadata_sequencing[metadata_sequencing$Libraries %in% i,]
}
metadata.df <- as.data.frame(matrix(nrow = 0, ncol = length(metadata[[1]])))
for (i in seq_along(metadata)) {
  meta.df <- as.data.frame(metadata[[i]])
  metadata.df <- rbind(metadata.df, meta.df)
}
remove(meta.df)
metadata.df = metadata.df[, c(2,6)]
orig.ident = as.data.frame(scrna_mimic$orig.ident)
colnames(orig.ident) = "sequencing.ident"
mimic.id = metadata.df$`Sample Labels`
for (i in seq_along(filenames)) {
  orig.ident[, "sequencing.ident"] <- gsub(filenames[i], mimic.id[i], orig.ident[, "sequencing.ident"])
}
colnames(orig.ident) = "mimic.id"
scrna_mimic = AddMetaData(scrna_mimic, metadata = orig.ident, col.name = "mimic.id")
scrna_mimic_merge = AddMetaData(scrna_mimic_merge, metadata = orig.ident, col.name = "mimic.id")
DimPlot(scrna_mimic_merge, group.by = "mimic.id", label = T, raster = F)
DimPlot(scrna_mimic, group.by = "mimic.id", label = T, raster = F)

```
```{r Cancer type}
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata_ID$`MIMIC ID` = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
mimic.id = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
metadata_ID = metadata_ID[, c(3,4)]
colnames(metadata_ID) = c("mimic.id", "Disease")
orig.ident = as.data.frame(cbind(scrna_mimic$mimic.id, scrna_mimic$orig.ident))
colnames(orig.ident) = c("mimic.id", "orig.ident")
disease.ident = merge(x = orig.ident, 
                   y = metadata_ID,
                   by = "mimic.id")
rownames(disease.ident) = rownames(orig.ident)
disease.ident[1:2] = NULL
scrna_mimic = AddMetaData(scrna_mimic, metadata = disease.ident, col.name = "Cancer.type")
DimPlot(scrna_mimic, group.by = "Cancer.type")

write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")

```
##Trajectory + pseudotime
```{r Trajectory and Pseudotime by monocle3}
scrna_mimic = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
cds_mimic = as.cell_data_set(scrna_mimic)
#partitions:
cds_mimic <- cluster_cells(cds_mimic, resolution=1e-3)
scrna_mimic_sub <- subset(as.Seurat(cds_mimic, assay = NULL), monocle3_partitions == 1)
scrna_mimic_sub = subset(scrna_mimic[,colnames(scrna_mimic)%in%colnames(scrna_mimic_sub)])
cds_mimic = as.cell_data_set(scrna_mimic_sub)
list_cluster <- scrna_mimic_sub@meta.data[["SCT_snn_res.0.8"]] ###IMPORTANT: Change to cell type when annotation is done
names(list_cluster) <- scrna_mimic_sub@assays[["SCT"]]@data@Dimnames[[2]]
cds_mimic@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_mimic@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
reducedDim(cds_mimic, "UMAP", withDimnames=TRUE) =scrna_mimic_sub@reductions[["umap"]]@cell.embeddings
cds_mimic@reduce_dim_aux$gene_loadings <- scrna_mimic_sub@reductions[["pca"]]@feature.loadings
recreate.partition <- c(rep(1, length(cds_mimic@colData@rownames)))
names(recreate.partition) <- cds_mimic@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_mimic@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
cds_mimic <- learn_graph(cds_mimic, 
                         use_partition = T, 
                         verbose = T, 
                         learn_graph_control=list(ncenter=1000),
                         close_loop = F)
plot_cells(cds_mimic,
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,group_label_size = 2,
           label_leaves=FALSE,
           label_branch_points=FALSE) #+ 
  #facet_wrap(~Cancer.type, nrow = 2)

cds_mimic <- order_cells(cds_mimic)
#save_monocle_objects(cds_mimic, directory_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Trajectory/cds_mimic")
plot_cells(cds_mimic,
           color_cells_by = "pseudotime", 
           cell_size = 1.25,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 1.5)


```
##Integration with fetal development dataset
```{r Integration Fetal}
scrna_mimic=read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]

raw_reads = list()
seurat_objects=list()
scrna_matrix = list()
matrix_list = list()
Assay = list()

Idents(scrna_kaesmann_subsampled) = "batch"
for (i in unique(scrna_kaesmann_subsampled$batch)) {
 raw_reads[[i]] = subset(scrna_kaesmann_subsampled, ident = i)
 genes=raw_reads[[i]]@assays[["RNA"]]@meta.features$feature_name
 metafeatures = raw_reads[[i]]@assays[["RNA"]]@meta.features
 metafeatures$ENSG = rownames(metafeatures)
 scrna_matrix[[i]] = raw_reads[[i]]@assays[["RNA"]]@counts
 rownames(raw_reads[[i]]@assays[["RNA"]]@counts) = ifelse(scrna_matrix[[i]]@Dimnames[[1]] %in%  metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
}
Idents(scrna_mimic) = "orig.ident"
for (i in unique(scrna_mimic$orig.ident)) {
 raw_reads[[i]] = subset(scrna_mimic, ident = i)
}
IDs = paste0(c(unique(as.character(scrna_kaesmann_subsampled$batch)), c(unique(as.character(scrna_mimic$orig.ident)))))
for (i in IDs) {
  matrix_list[[i]]= raw_reads[[i]]@assays$RNA$counts
  seurat_objects[[i]] = CreateSeuratObject(counts = matrix_list[[i]], project = i, min.cells = 3, min.features = 150) 
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay[[i]] = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[[i]][rb.genes,])/Matrix::colSums(Assay[[i]])*100
seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
}
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
scrna_merged = merge(seurat_objects_first, seurat_objects, merge.data = TRUE, project = "Integration")
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
scrna_merged= CellCycleScoring(scrna_merged, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_merged<-RunPCA (scrna_merged, verbose = FALSE)
scrna_integrated <- IntegrateLayers(object = scrna_merged, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_integrated <- FindNeighbors(object = scrna_integrated, reduction = "harmony", dims = 1:30)
scrna_integrated = RunUMAP(scrna_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_integrated <- FindClusters(scrna_integrated, resolution = i)

Idents(scrna_integrated) = "orig.ident"
IDs[stringr::str_detect(IDs, pattern = "SN")] = "Fetal"
IDs[stringr::str_detect(IDs, pattern = "LX")] = "All"
meta = as.data.frame(scrna_integrated$orig.ident)
meta$new.idents = ifelse(stringr::str_detect(meta$`scrna_integrated$orig.ident`, pattern = "SN"), "Fetal", "All")
meta$`scrna_integrated$orig.ident` = NULL
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "new.idents")
meta = as.data.frame(scrna_kaesmann_subsampled$precisest_label)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "precisest.label")
meta = as.data.frame(scrna_kaesmann_subsampled$dev_state)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "dev.state")
write_rds(scrna_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/Integration_fetal_all/scrna_integrated.rds")
```
```{r InferCNV_xCelltype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/Integration_fetal_all/scrna_integrated.rds")
Idents(scrna_integrated) = "anno_stage1"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "new.idents")
scrna_integrated_split$Fetal$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$Fetal), "_Fetal"))
DimPlot(scrna_integrated_split$Fetal, group.by = "anno_infercnv")
scrna_integrated_split$MB$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$MB), "_All"))
DimPlot(scrna_integrated_split$MB, group.by = "anno_infercnv")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(barcodes[[i]]) = "cell.type"
}
barcodes = Join(data.list = barcodes, by = "cell.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(annotation[[i]]) = "cell.type"
}
annotation = Join(data.list = annotation, by = "cell.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names=unique(as.character(scrna_integrated_split$Fetal$anno_infercnv)))

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCelltype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCelltype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")
```
```{r InferCNV_xCancertype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
Idents(scrna_integrated) = "cancer.type"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "cancer.type")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(barcodes[[i]]) = "cancer.type"
}
barcodes = Join(data.list = barcodes, by = "cancer.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(annotation[[i]]) = "cancer.type"
}
annotation = Join(data.list = annotation, by = "cancer.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names="Fetal")

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")

```

#Medulloblastoma
##Preprocess_MB
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r libraries, include=FALSE}
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(readr)
library(SingleR)
library(pheatmap)
library(dittoSeq)
library(SeuratWrappers)
library(monocle3)
library(DoubletFinder)
library(magrittr)
library(tradeSeq)
library(DropletUtils)
library(enrichR)
library(infercnv)
library(tidyr)
library(str2str)
library(NGCHM)
library(infercnvNGCHM)
library(futile.logger)
```
```{r MIMIC input files, include=FALSE}
base_directory = getwd() #Make sure you have the file you want to analyze in this folder: Name the folder with sequencing facility ID and put filtered_feature_bc_matrix in it, the script will fetch the data that it needs
filenames = c(paste0("LX", 430:432), paste0("LX", 434:436), paste0("LX", 438), paste0("LX", 440:441), paste0("LX", 500:504),paste0("LX", 507))
raw_reads = list()
seurat_objects=list()


for (i in filenames) {
  # Create the file name for each sample
  file_dir = paste0(base_directory,"/10x_runs/", i, "/filtered_feature_bc_matrix")
  # Read the sample data from the file using Read10x
  sample_data = Read10X(file_dir)
  raw_reads[[i]] = sample_data
  seurat_objects[[i]] = CreateSeuratObject(counts = raw_reads[[i]], project = i, min.cells = 3, min.features = 150)
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[rb.genes,])/Matrix::colSums(Assay)*100
  seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
  seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
  }
rm(sample_data, raw_reads) #Housekeeping

```
```{r QC, include=FALSE}
for (i in filenames) {
  min.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) - 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) + 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nCount.thr = median(seurat_objects[[i]]$nCount_RNA) + 3*mad(seurat_objects[[i]]$nCount_RNA)
  QC = seurat_objects[[i]]@meta.data
  QC$nCount_RNA = QC$nCount_RNA
  QC$nFeature_RNA = QC$nFeature_RNA
  QC$Dropouts = ifelse(QC$nCount_RNA<max.nCount.thr&
                         QC$nFeature_RNA>min.nFeature.thr&
                         QC$nFeature_RNA<max.nFeature.thr, FALSE, TRUE)
  p = QC %>%
    arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, colour=percent.mt, shape = Dropouts)) + 
    geom_point() + 
    scale_shape_manual(values = c(16,17), 
                       labels = c(paste0("Retained (",sum(QC$Dropouts == FALSE), " cells)"),
                                  paste0("Dropouts (",sum(QC$Dropouts == TRUE), " cells)"))) +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle(paste0("QC metrics ", unique(QC$orig.ident)))+
    geom_vline(xintercept = max.nCount.thr)+
    geom_hline(yintercept = c(min.nFeature.thr,max.nFeature.thr))+
    ylim(min.nFeature.thr-1,max(QC$nFeature_RNA))
  print(p)
}
rm(p,QC)
qc_seurat_objects <- list()
for (i in filenames) {
  min.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) - 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) + 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nCount.thr = median(seurat_objects[[i]]$nCount_RNA) + 3*mad(seurat_objects[[i]]$nCount_RNA)
  qc_seurat_objects[[i]] = subset(seurat_objects[[i]], subset = nFeature_RNA>min.nFeature.thr & 
                                 nFeature_RNA<max.nFeature.thr & 
                                 nCount_RNA<max.nCount.thr & 
                                 percent.mt<5 &
                                 percent.ribo<5)
}
rm(seurat_objects)


```
##Strategies_MB
```{r Data Merging -Important in Seurat v5- SCT Normalization, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects =  merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
seurat_objects= CellCycleScoring(seurat_objects, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_mimic_merge<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic_merge <- FindNeighbors(object = scrna_mimic_merge, reduction = "pca", dims = 1:30)
scrna_mimic_merge = RunUMAP(scrna_mimic_merge, reduction = "pca", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_mimic_merge <- FindClusters(scrna_mimic_merge, resolution = i)
DimPlot(scrna_mimic_merge, label = T, group.by = "orig.ident")
```
```{r Harmony Integration}
scrna_mimic<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = scrna_mimic, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
DimPlot(scrna_mimic, group.by = "orig.ident")
```

```{r Data Merging -Important in Seurat v5- SCT Normalization + CCA Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2")
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
write_rds(seurat_objects, "scrna_SCT_Normalized.rds")
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:30, reduction.name = "umap")
#write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_SCT_CCA/scrna_mb.rds")
```
```{r Data Merging -Important in Seurat v5- Log Normalization + Harmony Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects))#, vars.to.regress = c("percent.mt", "percent.ribo"))
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap")
#write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_Harmony/scrna_mb.rds")
DimPlot(scrna_mimic, split.by = "orig.ident")
```
```{r Data Merging -Important in Seurat v5- Log Normalization + CCA Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects), vars.to.regress = c("percent.mt", "percent.ribo"))
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:50)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:50, reduction.name = "umap")
#write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_CCA/scrna_mb.rds")
DimPlot(scrna_mimic, group.by = "orig.ident")
```
##Metadata curation_MB
```{r MimicID}
metadata_sequencing = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 1)
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata = list()
for (i in filenames) {
  metadata[[i]]= metadata_sequencing[metadata_sequencing$Libraries %in% i,]
}
metadata.df <- as.data.frame(matrix(nrow = 0, ncol = length(metadata[[1]])))
for (i in seq_along(metadata)) {
  meta.df <- as.data.frame(metadata[[i]])
  metadata.df <- rbind(metadata.df, meta.df)
}
remove(meta.df)
metadata.df = metadata.df[, c(2,6)]
orig.ident = as.data.frame(scrna_mimic$orig.ident)
colnames(orig.ident) = "sequencing.ident"
mimic.id = metadata.df$`Sample Labels`
for (i in seq_along(filenames)) {
  orig.ident[, "sequencing.ident"] <- gsub(filenames[i], mimic.id[i], orig.ident[, "sequencing.ident"])
}
colnames(orig.ident) = "mimic.id"
scrna_mimic = AddMetaData(scrna_mimic, metadata = orig.ident, col.name = "mimic.id")
scrna_mimic_merge = AddMetaData(scrna_mimic_merge, metadata = orig.ident, col.name = "mimic.id")
DimPlot(scrna_mimic_merge, group.by = "mimic.id", label = T, raster = F)
DimPlot(scrna_mimic, group.by = "mimic.id", label = T, raster = F)

```
```{r Cancer type}
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata_ID$`MIMIC ID` = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
mimic.id = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
metadata_ID = metadata_ID[, c(3,4)]
colnames(metadata_ID) = c("mimic.id", "Disease")
orig.ident = as.data.frame(cbind(scrna_mimic$mimic.id, scrna_mimic$orig.ident))
colnames(orig.ident) = c("mimic.id", "orig.ident")
disease.ident = merge(x = orig.ident, 
                   y = metadata_ID,
                   by = "mimic.id")
rownames(disease.ident) = rownames(orig.ident)
disease.ident[1:2] = NULL
scrna_mimic = AddMetaData(scrna_mimic, metadata = disease.ident, col.name = "Cancer.type")
DimPlot(scrna_mimic, group.by = "Cancer.type")

write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/scrna_mimic_mb.rds")

```
```{r Cell-type annotation through SingleR}
scrna_mimic = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/scrna_mimic_mb.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]
genes=scrna_kaesmann_subsampled@assays[["RNA"]]@meta.features$feature_name
metafeatures = scrna_kaesmann_subsampled@assays[["RNA"]]@meta.features
metafeatures$ENSG = rownames(metafeatures)
scrna_matrix = scrna_kaesmann_subsampled@assays[["RNA"]]@counts
rownames(scrna_kaesmann_subsampled@assays[["RNA"]]@counts) = ifelse(scrna_matrix@Dimnames[[1]] %in% metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
filenames.counts = paste0("counts.", filenames)
DefaultAssay(scrna_mimic) = "RNA"
scrna_mimic = JoinLayers(scrna_mimic)
annotations = SingleR(test=scrna_mimic@assays$RNA$counts,
                             ref=scrna_kaesmann_subsampled@assays$RNA$counts, labels=scrna_kaesmann_subsampled$precisest_label,
                             de.method="wilcox", de.n = 10)
plotScoreHeatmap(annotations)
transfer.anno = as.data.frame(annotations$labels, row.names = rownames(annotations))
transfer.anno$`annotations$labels` = as.factor(transfer.anno$`annotations$labels`)
scrna_mimic <- AddMetaData(scrna_mimic, transfer.anno, col.name = "precisestlabel")
DimPlot(scrna_mimic, reduction="umap", group.by = "precisestlabel", label = T) 
write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/scrna_mimic_mb.rds")
write_rds(annotations, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Annotations/annotations_analysis_fromSingleR.rds")
Idents(scrna_mimic) = "precisestlabel"
for (i in unique(Idents(scrna_mimic))) {
  scrna_mimic = AddMetaData(scrna_mimic, metadata = annotations$scores[,i], col.name = i)
}


```
```{r Doublet, eval=FALSE, include=FALSE}
remove.packages("Matrix")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
scrna_mimic = JoinLayers(scrna_mimic, assay = 'RNA', layers = 'counts')
sweep.res.mimic <- paramSweep(scrna_mimic, PCs = 1:30, sct = T)
sweep.stats <- summarizeSweep(sweep.res.mimic, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
annotations <- Idents(scrna_mimic)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi = 3200
print(paste0("Expected number of doublets: ", nExp_poi))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scrna_mimic <- doubletFinder(scrna_mimic, PCs = 1:30, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
doublets <- as.data.frame(cbind(colnames(scrna_mimic), scrna_mimic@meta.data[,grepl(paste0("pANN_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(scrna_mimic@meta.data))], scrna_mimic@meta.data[,grepl(paste0("DF.classifications_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(scrna_mimic@meta.data))]))
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)
rownames(doublets) = doublets$Barcode
doublets_meta = doublets
doublets_meta$DoubletFinder_score = NULL
doublets_meta$Barcode = NULL
scrna_mimic = AddMetaData(scrna_mimic, doublets_meta, "doublets")
DimPlot(scrna_mimic, split.by = "doublets")
```
##Trajectory + pseudotime_MB
```{r Trajectory and Pseudotime by monocle3}
scrna_mimic = read_rds("MB/scrna_mimic_mb.rds")
cds_mimic = as.cell_data_set(scrna_mimic)
#partitions:
cds_mimic <- cluster_cells(cds_mimic, resolution=1e-3)
scrna_mimic_sub <- subset(as.Seurat(cds_mimic, assay = NULL), monocle3_partitions == 1)
scrna_mimic_sub = subset(scrna_mimic[,colnames(scrna_mimic)%in%colnames(scrna_mimic_sub)])
cds_mimic = as.cell_data_set(scrna_mimic_sub)
list_cluster <- scrna_mimic_sub@meta.data[["SCT_snn_res.0.8"]] ###IMPORTANT: Change to cell type when annotation is done
names(list_cluster) <- scrna_mimic_sub@assays[["SCT"]]@data@Dimnames[[2]]
cds_mimic@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_mimic@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
reducedDim(cds_mimic, "UMAP", withDimnames=TRUE) =scrna_mimic_sub@reductions[["umap"]]@cell.embeddings
cds_mimic@reduce_dim_aux$gene_loadings <- scrna_mimic_sub@reductions[["pca"]]@feature.loadings
recreate.partition <- c(rep(1, length(cds_mimic@colData@rownames)))
names(recreate.partition) <- cds_mimic@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_mimic@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
cds_mimic <- learn_graph(cds_mimic, 
                         use_partition = T, 
                         verbose = T, 
                         learn_graph_control=list(ncenter=1000),
                         close_loop = F)
plot_cells(cds_mimic,
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,group_label_size = 2,
           label_leaves=FALSE,
           label_branch_points=FALSE) #+ 
  #facet_wrap(~Cancer.type, nrow = 2)

cds_mimic <- order_cells(cds_mimic)
#save_monocle_objects(cds_mimic, directory_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Trajectory/cds_mimic")
plot_cells(cds_mimic,
           color_cells_by = "pseudotime", 
           cell_size = 1.25,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 1.5)


```
##Integration with fetal development dataset_MB
```{r Integration Fetal and Medulloblastoma}
scrna_mimic=read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/scrna_mimic_mb.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]

raw_reads = list()
seurat_objects=list()
scrna_matrix = list()
matrix_list = list()
Assay = list()

Idents(scrna_kaesmann_subsampled) = "batch"
for (i in unique(scrna_kaesmann_subsampled$batch)) {
 raw_reads[[i]] = subset(scrna_kaesmann_subsampled, ident = i)
 genes=raw_reads[[i]]@assays[["RNA"]]@meta.features$feature_name
 metafeatures = raw_reads[[i]]@assays[["RNA"]]@meta.features
 metafeatures$ENSG = rownames(metafeatures)
 scrna_matrix[[i]] = raw_reads[[i]]@assays[["RNA"]]@counts
 rownames(raw_reads[[i]]@assays[["RNA"]]@counts) = ifelse(scrna_matrix[[i]]@Dimnames[[1]] %in%  metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
}
Idents(scrna_mimic) = "orig.ident"
for (i in unique(scrna_mimic$orig.ident)) {
 raw_reads[[i]] = subset(scrna_mimic, ident = i)
}
IDs = paste0(c(unique(as.character(scrna_kaesmann_subsampled$batch)), c(unique(as.character(scrna_mimic$orig.ident)))))
for (i in IDs) {
  matrix_list[[i]]= raw_reads[[i]]@assays$RNA$counts
  seurat_objects[[i]] = CreateSeuratObject(counts = matrix_list[[i]], project = i, min.cells = 3, min.features = 150) 
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay[[i]] = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[[i]][rb.genes,])/Matrix::colSums(Assay[[i]])*100
seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
}
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
scrna_merged = merge(seurat_objects_first, seurat_objects, merge.data = TRUE, project = "Integration")
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
scrna_merged= CellCycleScoring(scrna_merged, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_merged<-RunPCA (scrna_merged, verbose = FALSE)
scrna_integrated <- IntegrateLayers(object = scrna_merged, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_integrated <- FindNeighbors(object = scrna_integrated, reduction = "harmony", dims = 1:30)
scrna_integrated = RunUMAP(scrna_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_integrated <- FindClusters(scrna_integrated, resolution = i)

Idents(scrna_integrated) = "orig.ident"
IDs[stringr::str_detect(IDs, pattern = "SN")] = "Fetal"
IDs[stringr::str_detect(IDs, pattern = "LX")] = "MB"
meta = as.data.frame(scrna_integrated$orig.ident)
meta$new.idents = ifelse(stringr::str_detect(meta$`scrna_integrated$orig.ident`, pattern = "SN"), "Fetal", "MB")
meta$`scrna_integrated$orig.ident` = NULL
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "new.idents")
meta = as.data.frame(scrna_kaesmann_subsampled$precisest_label)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "precisest.label")
meta = as.data.frame(scrna_kaesmann_subsampled$dev_state)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "dev.state")
write_rds(scrna_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
```
```{r InferCNV_xCelltype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
Idents(scrna_integrated) = "anno_stage1"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "new.idents")
scrna_integrated_split$Fetal$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$Fetal), "_Fetal"))
DimPlot(scrna_integrated_split$Fetal, group.by = "anno_infercnv")
scrna_integrated_split$MB$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$MB), "_MB"))
DimPlot(scrna_integrated_split$MB, group.by = "anno_infercnv")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(barcodes[[i]]) = "cell.type"
}
barcodes = Join(data.list = barcodes, by = "cell.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(annotation[[i]]) = "cell.type"
}
annotation = Join(data.list = annotation, by = "cell.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names=unique(as.character(scrna_integrated_split$Fetal$anno_infercnv)))

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCelltype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCelltype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")
```
```{r InferCNV_xCancertype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
Idents(scrna_integrated) = "cancer.type"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "cancer.type")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(barcodes[[i]]) = "cancer.type"
}
barcodes = Join(data.list = barcodes, by = "cancer.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(annotation[[i]]) = "cancer.type"
}
annotation = Join(data.list = annotation, by = "cancer.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names="Fetal")

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")

```

#Ependymoma
##Preprocess_EP
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
```{r libraries, include=FALSE}
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(readr)
library(SingleR)
library(pheatmap)
library(dittoSeq)
library(SeuratWrappers)
library(monocle3)
library(DoubletFinder)
library(magrittr)
library(tradeSeq)
library(DropletUtils)
library(enrichR)
library(infercnv)
library(tidyr)
library(str2str)
library(NGCHM)
library(infercnvNGCHM)
library(futile.logger)
```
```{r MIMIC input files, include=FALSE}
base_directory = getwd() #Make sure you have the file you want to analyze in this folder: Name the folder with sequencing facility ID and put filtered_feature_bc_matrix in it, the script will fetch the data that it needs
filenames = c(paste0("LX", 433), paste0("LX", 437), paste0("LX", 439), paste0("LX", 505:506))
raw_reads = list()
seurat_objects=list()


for (i in filenames) {
  # Create the file name for each sample
  file_dir = paste0(base_directory,"/10x_runs/", i, "/filtered_feature_bc_matrix")
  # Read the sample data from the file using Read10x
  sample_data = Read10X(file_dir)
  raw_reads[[i]] = sample_data
  seurat_objects[[i]] = CreateSeuratObject(counts = raw_reads[[i]], project = i, min.cells = 3, min.features = 150)
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[rb.genes,])/Matrix::colSums(Assay)*100
  seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
  seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
  }
rm(sample_data, raw_reads) #Housekeeping

```
```{r QC, include=FALSE}
for (i in filenames) {
  min.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) - 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) + 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nCount.thr = median(seurat_objects[[i]]$nCount_RNA) + 3*mad(seurat_objects[[i]]$nCount_RNA)
  QC = seurat_objects[[i]]@meta.data
  QC$nCount_RNA = QC$nCount_RNA
  QC$nFeature_RNA = QC$nFeature_RNA
  QC$Dropouts = ifelse(QC$nCount_RNA<max.nCount.thr&
                         QC$nFeature_RNA>min.nFeature.thr&
                         QC$nFeature_RNA<max.nFeature.thr, FALSE, TRUE)
  p = QC %>%
    arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, colour=percent.mt, shape = Dropouts)) + 
    geom_point() + 
    scale_shape_manual(values = c(16,17), 
                       labels = c(paste0("Retained (",sum(QC$Dropouts == FALSE), " cells)"),
                                  paste0("Dropouts (",sum(QC$Dropouts == TRUE), " cells)"))) +
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
    ggtitle(paste0("QC metrics ", unique(QC$orig.ident)))+
    geom_vline(xintercept = max.nCount.thr)+
    geom_hline(yintercept = c(min.nFeature.thr,max.nFeature.thr))+
    ylim(min.nFeature.thr-1,max(QC$nFeature_RNA))
  print(p)
}
rm(p,QC)
qc_seurat_objects <- list()
for (i in filenames) {
  min.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) - 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nFeature.thr = median(seurat_objects[[i]]$nFeature_RNA) + 3*mad(seurat_objects[[i]]$nFeature_RNA)
  max.nCount.thr = median(seurat_objects[[i]]$nCount_RNA) + 3*mad(seurat_objects[[i]]$nCount_RNA)
  qc_seurat_objects[[i]] = subset(seurat_objects[[i]], subset = nFeature_RNA>min.nFeature.thr & 
                                 nFeature_RNA<max.nFeature.thr & 
                                 nCount_RNA<max.nCount.thr & 
                                 percent.mt<5 &
                                 percent.ribo<5)
}
rm(seurat_objects)


```
##Strategies_EP
```{r Data Merging -Important in Seurat v5- SCT Normalization, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects =  merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
seurat_objects= CellCycleScoring(seurat_objects, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_mimic_merge<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic_merge <- FindNeighbors(object = scrna_mimic_merge, reduction = "pca", dims = 1:30)
scrna_mimic_merge = RunUMAP(scrna_mimic_merge, reduction = "pca", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_mimic_merge <- FindClusters(scrna_mimic_merge, resolution = i)
DimPlot(scrna_mimic_merge, label = T, group.by = "orig.ident")
```
```{r Harmony Integration}
scrna_mimic<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = scrna_mimic, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
DimPlot(scrna_mimic, group.by = "orig.ident")
```

```{r Data Merging -Important in Seurat v5- SCT Normalization + CCA Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2")
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:30, reduction.name = "umap")
#write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_SCT_CCA/scrna.rds")
```
```{r Data Merging -Important in Seurat v5- Log Normalization + Harmony Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects))#, vars.to.regress = c("percent.mt", "percent.ribo"))
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "harmony", dims = 1:30)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "harmony", dims = 1:30, reduction.name = "umap")
#write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_Harmony/scrna.rds")
DimPlot(scrna_mimic, split.by = "orig.ident")
```
```{r Data Merging -Important in Seurat v5- Log Normalization + CCA Integration, eval=FALSE, include=FALSE}
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects = merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
seurat_objects = NormalizeData(seurat_objects, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_objects = FindVariableFeatures(seurat_objects, selection.method = "vst", nfeatures = 2000)
seurat_objects = ScaleData(seurat_objects, features = rownames(seurat_objects), vars.to.regress = c("percent.mt", "percent.ribo"))
seurat_objects<-RunPCA (seurat_objects, verbose = FALSE)
scrna_mimic <- IntegrateLayers(object = seurat_objects, method = CCAIntegration, orig.reduction = "pca", normalization.method = "LogNormalize")
scrna_mimic <- FindNeighbors(object = scrna_mimic, reduction = "integrated.dr", dims = 1:50)
i = seq(0.2, 1, by = 0.2)
scrna_mimic <- FindClusters(scrna_mimic, resolution = i)
scrna_mimic = RunUMAP(scrna_mimic, reduction = "integrated.dr", dims = 1:50, reduction.name = "umap")
#write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/scrna_Log_CCA/scrna.rds")
DimPlot(scrna_mimic, group.by = "orig.ident")
```
##Metadata curation_EP
```{r MimicID}
metadata_sequencing = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 1)
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata = list()
for (i in filenames) {
  metadata[[i]]= metadata_sequencing[metadata_sequencing$Libraries %in% i,]
}
metadata.df <- as.data.frame(matrix(nrow = 0, ncol = length(metadata[[1]])))
for (i in seq_along(metadata)) {
  meta.df <- as.data.frame(metadata[[i]])
  metadata.df <- rbind(metadata.df, meta.df)
}
remove(meta.df)
metadata.df = metadata.df[, c(2,6)]
orig.ident = as.data.frame(scrna_mimic$orig.ident)
colnames(orig.ident) = "sequencing.ident"
mimic.id = metadata.df$`Sample Labels`
for (i in seq_along(filenames)) {
  orig.ident[, "sequencing.ident"] <- gsub(filenames[i], mimic.id[i], orig.ident[, "sequencing.ident"])
}
colnames(orig.ident) = "mimic.id"
scrna_mimic = AddMetaData(scrna_mimic, metadata = orig.ident, col.name = "mimic.id")
scrna_mimic_merge = AddMetaData(scrna_mimic_merge, metadata = orig.ident, col.name = "mimic.id")
DimPlot(scrna_mimic_merge, group.by = "mimic.id", label = T, raster = F)
DimPlot(scrna_mimic, group.by = "mimic.id", label = T, raster = F)

```
```{r Cancer type}
metadata_ID = readxl::read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Metadata.xlsx", sheet = 2)
metadata_ID$`MIMIC ID` = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
mimic.id = stringr::str_replace(metadata_ID$`MIMIC ID`, pattern = "-", replacement = '')
metadata_ID = metadata_ID[, c(3,4)]
colnames(metadata_ID) = c("mimic.id", "Disease")
orig.ident = as.data.frame(cbind(scrna_mimic$mimic.id, scrna_mimic$orig.ident))
colnames(orig.ident) = c("mimic.id", "orig.ident")
disease.ident = merge(x = orig.ident, 
                   y = metadata_ID,
                   by = "mimic.id")
rownames(disease.ident) = rownames(orig.ident)
disease.ident[1:2] = NULL
scrna_mimic = AddMetaData(scrna_mimic, metadata = disease.ident, col.name = "Cancer.type")
DimPlot(scrna_mimic, split.by = "Cancer.type")

write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/scrna_mimic_ep.rds")

```
```{r Cell-type annotation through SingleR}
scrna_mimic = read_rds("EP/scrna_mimic_ep.rds")
scrna_kaesmann = readRDS(file = "~/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]
genes=scrna_kaesmann_subsampled@assays[["RNA"]]@meta.features$feature_name
metafeatures = scrna_kaesmann_subsampled@assays[["RNA"]]@meta.features
metafeatures$ENSG = rownames(metafeatures)
scrna_matrix = scrna_kaesmann_subsampled@assays[["RNA"]]@counts
rownames(scrna_kaesmann_subsampled@assays[["RNA"]]@counts) = ifelse(scrna_matrix@Dimnames[[1]] %in% metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
filenames.counts = paste0("counts.", filenames)
DefaultAssay(scrna_mimic) = "RNA"
scrna_mimic = JoinLayers(scrna_mimic)
annotations = SingleR(test=scrna_mimic@assays$RNA$counts,
                             ref=scrna_kaesmann_subsampled@assays$RNA$counts, labels=scrna_kaesmann_subsampled$precisest_label,
                             de.method="wilcox", de.n = 10)
plotScoreHeatmap(annotations)
transfer.anno = as.data.frame(annotations$labels, row.names = rownames(annotations))
transfer.anno$`annotations$labels` = as.factor(transfer.anno$`annotations$labels`)
scrna_mimic <- AddMetaData(scrna_mimic, transfer.anno, col.name = "precisestlabel")
DimPlot(scrna_mimic, reduction="umap", group.by = "precisestlabel", label = T) 
write_rds(scrna_mimic, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/scrna_mimic_ep.rds")
write_rds(annotations, "hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/Annotations/annotations_analysis_fromSingleR.rds")
Idents(scrna_mimic) = "precisestlabel"
for (i in unique(Idents(scrna_mimic))) {
  scrna_mimic = AddMetaData(scrna_mimic, metadata = annotations$scores[,i], col.name = i)
}


```
```{r Doublet, eval=FALSE, include=FALSE}
remove.packages("Matrix")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
scrna_mimic = JoinLayers(scrna_mimic, assay = 'RNA', layers = 'counts')
sweep.res.mimic <- paramSweep(scrna_mimic, PCs = 1:30, sct = T)
sweep.stats <- summarizeSweep(sweep.res.mimic, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
annotations <- Idents(scrna_mimic)
homotypic.prop <- modelHomotypic(annotations)
nExp_poi = 3200
print(paste0("Expected number of doublets: ", nExp_poi))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scrna_mimic <- doubletFinder(scrna_mimic, PCs = 1:30, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
doublets <- as.data.frame(cbind(colnames(scrna_mimic), scrna_mimic@meta.data[,grepl(paste0("pANN_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(scrna_mimic@meta.data))], scrna_mimic@meta.data[,grepl(paste0("DF.classifications_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(scrna_mimic@meta.data))]))
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)
rownames(doublets) = doublets$Barcode
doublets_meta = doublets
doublets_meta$DoubletFinder_score = NULL
doublets_meta$Barcode = NULL
scrna_mimic = AddMetaData(scrna_mimic, doublets_meta, "doublets")
DimPlot(scrna_mimic, split.by = "doublets")
```
##Trajectory + pseudotime_EP
```{r Trajectory and Pseudotime by monocle3}
scrna_mimic = read_rds("hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/scrna_mimic_ep.rds")
cds_mimic = as.cell_data_set(scrna_mimic)
#partitions:
cds_mimic <- cluster_cells(cds_mimic, resolution=1e-3)
scrna_mimic_sub <- subset(as.Seurat(cds_mimic, assay = NULL), monocle3_partitions == 1)
scrna_mimic_sub = subset(scrna_mimic[,colnames(scrna_mimic)%in%colnames(scrna_mimic_sub)])
cds_mimic = as.cell_data_set(scrna_mimic_sub)
list_cluster <- scrna_mimic_sub@meta.data[["SCT_snn_res.0.8"]] ###IMPORTANT: Change to cell type when annotation is done
names(list_cluster) <- scrna_mimic_sub@assays[["SCT"]]@data@Dimnames[[2]]
cds_mimic@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_mimic@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
reducedDim(cds_mimic, "UMAP", withDimnames=TRUE) =scrna_mimic_sub@reductions[["umap"]]@cell.embeddings
cds_mimic@reduce_dim_aux$gene_loadings <- scrna_mimic_sub@reductions[["pca"]]@feature.loadings
recreate.partition <- c(rep(1, length(cds_mimic@colData@rownames)))
names(recreate.partition) <- cds_mimic@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_mimic@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
cds_mimic <- learn_graph(cds_mimic, 
                         use_partition = T, 
                         verbose = T, 
                         learn_graph_control=list(ncenter=1000),
                         close_loop = F)
plot_cells(cds_mimic,
           color_cells_by = "cluster",
           label_groups_by_cluster=TRUE,group_label_size = 2,
           label_leaves=FALSE,
           label_branch_points=FALSE) #+ 
  #facet_wrap(~Cancer.type, nrow = 2)

cds_mimic <- order_cells(cds_mimic)
#save_monocle_objects(cds_mimic, directory_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Trajectory/cds_mimic")
plot_cells(cds_mimic,
           color_cells_by = "pseudotime", 
           cell_size = 1.25,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           trajectory_graph_color = "green",
           trajectory_graph_segment_size = 1.5)


```
##Integration with fetal development dataset_EP
```{r Integration Fetal and Ependymoma}
scrna_mimic=read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/scrna_mimic_ep.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]

raw_reads = list()
seurat_objects=list()
scrna_matrix = list()
matrix_list = list()
Assay = list()

Idents(scrna_kaesmann_subsampled) = "batch"
for (i in unique(scrna_kaesmann_subsampled$batch)) {
 raw_reads[[i]] = subset(scrna_kaesmann_subsampled, ident = i)
 genes=raw_reads[[i]]@assays[["RNA"]]@meta.features$feature_name
 metafeatures = raw_reads[[i]]@assays[["RNA"]]@meta.features
 metafeatures$ENSG = rownames(metafeatures)
 scrna_matrix[[i]] = raw_reads[[i]]@assays[["RNA"]]@counts
 rownames(raw_reads[[i]]@assays[["RNA"]]@counts) = ifelse(scrna_matrix[[i]]@Dimnames[[1]] %in%  metafeatures$ENSG, as.character(metafeatures$feature_name), NA)
}
Idents(scrna_mimic) = "orig.ident"
for (i in unique(scrna_mimic$orig.ident)) {
 raw_reads[[i]] = subset(scrna_mimic, ident = i)
}
IDs = paste0(c(unique(as.character(scrna_kaesmann_subsampled$batch)), c(unique(as.character(scrna_mimic$orig.ident)))))
for (i in IDs) {
  matrix_list[[i]]= raw_reads[[i]]@assays$RNA$counts
  seurat_objects[[i]] = CreateSeuratObject(counts = matrix_list[[i]], project = i, min.cells = 3, min.features = 150) 
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))]
  Assay[[i]] = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[[i]][rb.genes,])/Matrix::colSums(Assay[[i]])*100
seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-")
seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo")
}
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
scrna_merged = merge(seurat_objects_first, seurat_objects, merge.data = TRUE, project = "Integration")
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
scrna_merged= CellCycleScoring(scrna_merged, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
scrna_merged = SCTransform(scrna_merged, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
scrna_merged<-RunPCA (scrna_merged, verbose = FALSE)
scrna_integrated <- IntegrateLayers(object = scrna_merged, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
scrna_integrated <- FindNeighbors(object = scrna_integrated, reduction = "harmony", dims = 1:30)
scrna_integrated = RunUMAP(scrna_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_integrated <- FindClusters(scrna_integrated, resolution = i)

Idents(scrna_integrated) = "orig.ident"
IDs[stringr::str_detect(IDs, pattern = "SN")] = "Fetal"
IDs[stringr::str_detect(IDs, pattern = "LX")] = "EP"
meta = as.data.frame(scrna_integrated$orig.ident)
meta$new.idents = ifelse(stringr::str_detect(meta$`scrna_integrated$orig.ident`, pattern = "SN"), "Fetal", "EP")
meta$`scrna_integrated$orig.ident` = NULL
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "new.idents")
meta = as.data.frame(scrna_kaesmann_subsampled$precisest_label)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "precisest.label")
meta = as.data.frame(scrna_kaesmann_subsampled$dev_state)
scrna_integrated=AddMetaData(scrna_integrated, metadata = meta, col.name = "dev.state")
write_rds(scrna_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/Integration_fetal_ep/scrna_integrated.rds")
```
```{r InferCNV_xCelltype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/Integration_fetal_ep/scrna_integrated.rds")
Idents(scrna_integrated) = "anno_stage1"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "new.idents")
scrna_integrated_split$Fetal$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$Fetal), "_Fetal"))
DimPlot(scrna_integrated_split$Fetal, group.by = "anno_infercnv")
scrna_integrated_split$MB$anno_infercnv = as.factor(paste0(Idents(scrna_integrated_split$MB), "_EP"))
DimPlot(scrna_integrated_split$MB, group.by = "anno_infercnv")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(barcodes[[i]]) = "cell.type"
}
barcodes = Join(data.list = barcodes, by = "cell.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$anno_infercnv)
  colnames(annotation[[i]]) = "cell.type"
}
annotation = Join(data.list = annotation, by = "cell.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names=unique(as.character(scrna_integrated_split$Fetal$anno_infercnv)))

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCelltype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCelltype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")
```
```{r InferCNV_xCancertype}
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/Integration_fetal_ep/scrna_integrated.rds")
Idents(scrna_integrated) = "cancer.type"
counts_matrix = GetAssayData(scrna_integrated, layer="counts")
scrna_integrated_split = SplitObject(scrna_integrated, split.by = "cancer.type")
barcodes = list()
for (i in seq_along(scrna_integrated_split)) {
  barcodes[[i]]= as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(barcodes[[i]]) = "cancer.type"
}
barcodes = Join(data.list = barcodes, by = "cancer.type")
annotation = list()
for (i in seq_along(scrna_integrated_split)) {
  annotation[[i]] = as.data.frame(scrna_integrated_split[[i]]$cancer.type)
  colnames(annotation[[i]]) = "cancer.type"
}
annotation = Join(data.list = annotation, by = "cancer.type")
annotation = as.matrix(annotation)
infercnv_mimic = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names="Fetal")

infercnv_mimic = infercnv::run(infercnv_mimic,
                             cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCancertype", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna_integrated = infercnv::add_to_seurat(seurat_obj = scrna_integrated,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCancertype",
                                           top_n = 15)
write_rds(infercnv_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/EP/infercnv/CNV_output_xCancertype/scrna_integrated_infercnv.rds")

```


#Annotation
```{r Annotations integrated dataset}
scrna_mimic=read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/scrna_mimic_mb.rds")
scrna_kaesmann = readRDS(file = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum Development/Kaesmann/seurat.rds")
scrna_kaesmann@assays[["RNA"]]@counts = scrna_kaesmann@assays[["RNA"]]@data  
scrna_kaesmann_subsampled <- scrna_kaesmann[, sample(colnames(scrna_kaesmann), size =ncol(scrna_mimic), replace=F)]
scrna_integrated = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
Idents(scrna_integrated)= "new.idents"
IDs = paste0(c(unique(as.character(scrna_kaesmann_subsampled$batch)), c(unique(as.character(scrna_mimic$orig.ident)))))
IDs[stringr::str_detect(IDs, pattern = "SN")] = "Fetal"
IDs[stringr::str_detect(IDs, pattern = "LX")] = "MB"

scrna_integrated_split = list()
for (i in unique(IDs)) {
  scrna_integrated_split[[i]] = subset(scrna_integrated, idents = i)
  DefaultAssay(scrna_integrated_split[[i]]) = "RNA"
  scrna_integrated_split[[i]] = JoinLayers(scrna_integrated_split[[i]])
}

annotations.precisest = SingleR(test=scrna_integrated_split$MB@assays$RNA$counts,
                             ref=scrna_integrated_split$Fetal@assays$RNA$counts,
                             labels=scrna_integrated_split$Fetal$precisest.label,
                             de.method="wilcox", de.n = 10)
write_rds(annotations.precisest, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/Annotations/annotations_integrated_fetal_mb.precisest.rds")
transfer.anno = as.data.frame(annotations.precisest$labels, row.names = rownames(annotations.precisest))
transfer.anno$`annotations.precisest$labels` = as.factor(transfer.anno$`annotations.precisest$labels`)
scrna_integrated_split$MB <- AddMetaData(scrna_integrated_split$MB, transfer.anno, col.name = "precisest.label")
barcodes= as.data.frame(c(colnames(scrna_integrated_split$Fetal), colnames(scrna_integrated_split$MB)))
metadata = as.data.frame(c(as.character(scrna_integrated_split$Fetal$precisest.label), as.character(scrna_integrated_split$MB$precisest.label)), row.names = barcodes$`c(colnames(scrna_integrated_split$Fetal), colnames(scrna_integrated_split$MB))`)
scrna_integrated=AddMetaData(scrna_integrated, metadata, "precisest.label")
#tab <- table(Assigned=annotations.precisest$labels, Cluster=scrna_integrated_split$MB$SCT_snn_res.0.8)
#pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))

annotations.dev.state = SingleR(test=scrna_integrated_split$MB@assays$RNA$counts,
                             ref=scrna_integrated_split$Fetal@assays$RNA$counts,
                             labels=scrna_integrated_split$Fetal$dev.state,
                             de.method="wilcox", de.n = 10)
write_rds(annotations.dev.state, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/Annotations/annotations_integrated_fetal_mb_devstate.rds")
transfer.anno = as.data.frame(annotations.dev.state$labels, row.names = rownames(annotations.dev.state))
transfer.anno$`annotations.dev.state$labels` = as.factor(transfer.anno$`annotations.dev.state$labels`)
scrna_integrated_split$MB <- AddMetaData(scrna_integrated_split$MB, transfer.anno, col.name = "dev.state")
barcodes= as.data.frame(c(colnames(scrna_integrated_split$Fetal), colnames(scrna_integrated_split$MB)))
metadata = as.data.frame(c(as.character(scrna_integrated_split$Fetal$dev.state), as.character(scrna_integrated_split$MB$dev.state)), row.names = barcodes$`c(colnames(scrna_integrated_split$Fetal), colnames(scrna_integrated_split$MB))`)
scrna_integrated=AddMetaData(scrna_integrated, metadata, "dev.state")
#tab <- table(Assigned=annotations.dev.state$labels, Cluster=scrna_integrated_split$MB$SCT_snn_res.0.8)
#pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
cancertype = as.data.frame(scrna_mimic$Cancer.type)
scrna_integrated = AddMetaData(scrna_integrated, cancertype, "cancer.type")
scrna_integrated$cancer.type = replace_na(scrna_integrated$cancer.type, "Fetal")
write_rds(scrna_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")


Meta = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/Annotations/Annotations_kaesmann.csv", sep = ",")
Meta$Gene_name = toupper(Meta$Gene_name)
Meta_hu_celltype = list()
for (i in unique(Meta$Cell.state)) {
  Meta_hu_celltype[[i]] = Meta[Meta$Cell.state == i,]
}
for (i in unique(Meta$Cell.state)) {
  scrna_integrated = AddModuleScore(scrna_integrated, list(Meta_hu_celltype[[i]]$Gene_name), name = paste0(i, "_score"))
}
for (i in unique(Meta$Cell.state)) {
  p = FeaturePlot(scrna_integrated, features = paste0(i, "_score1"), order=T, split.by = "new.idents")
  print(p)
}


Idents(scrna_integrated) = "SCT_snn_res.1"
scrna_integrated= PrepSCTFindMarkers(scrna_integrated)
DEG = FindAllMarkers(scrna_integrated, only.pos = T, min.diff.pct = 0.1)
DEG = DEG[DEG$p_val_adj<0.05,]
Idents(scrna_integrated) = "SCT_snn_res.1"
scrna_integrated = RenameIdents(scrna_integrated, c("0" = "GC/UBC",
                                                    "1" = "Purkinje",
                                                    "2" = "Progenitor_1",
                                                    "3" = "GC_diff_1",
                                                    "4" = "4",
                                                    "5" = "Progenitor_2",
                                                    "6" = "GC_diff_2",
                                                    "7" = "7",
                                                    "8" = "UBC",
                                                    "9" = "GC_diff_3",
                                                    "10" = "10",
                                                    "11" = "UBC",
                                                    "12" = "Immune",
                                                    "13" = "VZ_Neuroblast",
                                                    "14" = "GC_1",
                                                    "15" = "GC_diff_4",
                                                    "16" = "GC_diff_5",
                                                    "17" = "VZ_Neuroblast",
                                                    "18" = "GC_diff_6",
                                                    "19" = "Glut_DN",
                                                    "20" = "GC_diff_7",
                                                    "21" = "GABA_DN",
                                                    "22" = "22",
                                                    "23" = "EC_1",
                                                    "24" = "Oligodendrocyte",
                                                    "25" = "GC_2",
                                                    "26" = "Meningeal",
                                                    "27" = "GC_3",
                                                    "28" = "28",
                                                    "29" = "Oligodendrocyte",
                                                    "30" = "Interneuron_diff",
                                                    "31" = "Progenitor_2",
                                                    "32" = "Interneuron",
                                                    "33" = "33",
                                                    "34" = "EC_2",
                                                    "35" = "UBC",
                                                    "36" = "Progenitor_oligo",
                                                    "37" = "Bergmann",
                                                    "38" = "38", #Schwann Cells? 
                                                    "39" = "GABA_DN",
                                                    "40" = "Interneuron",
                                                    "41" = "41", #Podocytes?
                                                    "42" = "Astrocyte",
                                                    "43" = "GC_4"))


scrna_integrated$anno_stage1 = scrna_integrated@active.ident
write_rds(scrna_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/MB/Integration_fetal_mb/scrna_integrated.rds")
```



