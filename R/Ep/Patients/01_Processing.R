#Load libs
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(readr)
library(SingleR)
library(pheatmap)
library(dittoSeq)
library(SeuratWrappers)
library(DoubletFinder)
library(magrittr)
library(tradeSeq)
library(DropletUtils)
library(enrichR)
library(tidyr)
library(str2str)
library(NGCHM)
library(futile.logger)
library(scater)
#Setup
base_directory = getwd() #Make sure you have the file you want to analyze in this folder: Name the folder with sequencing facility ID and put filtered_feature_bc_matrix in it, the script will fetch the data that it needs
filenames = c("LX442", "LX498", "LX499") #Load the correct file labels
raw_reads = list() #Prepare raw reads list for for loop
seurat_objects=list() #Prepare seurat objects list for loop


#Objects setup
for (i in filenames) {
  # Create the file name for each sample
  file_dir = paste0(base_directory,"/fvalzano/Rstudio_Test1/Ependymoma_Patients/10x_runs/", i, "/filtered_feature_bc_matrix")
  # Read the sample data from the file using Read10x
  sample_data = Read10X(file_dir)
  raw_reads[[i]] = sample_data
  seurat_objects[[i]] = CreateSeuratObject(counts = raw_reads[[i]], project = i, min.cells = 3, min.features = 150)
  rb.genes = rownames(seurat_objects[[i]])[grep("^RP[SL]",rownames(seurat_objects[[i]]))] #Find ribosomial genes
  Assay = GetAssayData(seurat_objects[[i]])
  percent.ribo = colSums(Assay[rb.genes,])/Matrix::colSums(Assay)*100 #Calculate ribosomial protein gene content
  seurat_objects[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_objects[[i]], pattern = "^MT-") #Calculate mithochondrial gene content
  seurat_objects[[i]] = AddMetaData(seurat_objects[[i]], percent.ribo, col.name = "percent.ribo") #Manually add % of ribosomial protein genes
  }
rm(sample_data, raw_reads) #Housekeeping


#QC
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/QC.pdf", width=5, height=5)
for (i in filenames) {
  #calculation of threshold is performed either by either IsOutlier() function or via expliciting the formula - results are the same
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
   qc_seurat_objects[[i]] = subset(seurat_objects[[i]], subset = percent.mt <5 &
                                  percent.ribo<5)
}
rm(seurat_objects)

#SCT_Normalization
#In Seurat v5 first merge the objects and then integrate them
seurat_objects=qc_seurat_objects
seurat_objects_first = seurat_objects[[1]]
seurat_objects[[1]] = NULL
seurat_objects =  merge(x = seurat_objects_first, y= c(seurat_objects), merge.data = TRUE, project = "MIMIC") 
#Perform SCT normalization and then calculate cellcycle scoring, after that, rerun SCT regressing out the cell cycle Phase
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo"), vst.flavor = "v2", assay = "RNA")
seurat_objects= CellCycleScoring(seurat_objects, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, assay = 'SCT')
seurat_objects = SCTransform(seurat_objects, vars.to.regress = c("percent.mt", "percent.ribo", 'S.Score', 'G2M.Score'), vst.flavor = "v2", assay = "RNA")
#Perform dimensionality reduction on unintegrated object
scrna_ep<-RunPCA (seurat_objects, verbose = FALSE)
scrna_ep <- FindNeighbors(object = scrna_ep, reduction = "pca", dims = 1:30)
scrna_ep = RunUMAP(scrna_ep, reduction = "pca", dims = 1:30, reduction.name = "umap")
i = seq(0.2, 1, by = 0.2)
scrna_ep <- FindClusters(scrna_ep, resolution = i)
#Perform integration/batch correction
scrna_ep <- IntegrateLayers(object = scrna_ep, method = HarmonyIntegration, orig.reduction = "pca", normalization.method = "SCT")
#Perform dimensionality reduction on integrated object
scrna_ep <- FindNeighbors(object = scrna_ep, reduction = "harmony", dims = 1:30)
scrna_ep = RunUMAP(scrna_ep, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
i = seq(0.2, 1, by = 0.2)
scrna_ep <- FindClusters(scrna_ep, resolution = i)
DimPlot(scrna_ep, group.by = "orig.ident", reduction="umap_harmony")
write_rds(scrna_ep, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_Patients/scrna_ependymoma_all.rds")
