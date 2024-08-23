library(Seurat)
library(readr)
library(dplyr)

samples = c("D1_organoid", "D2_organoid", "D3_organoid")
files_10x = list()
Atamian_organoids = list()
#Set up the objects
for (s in samples) {
  files_10x[[s]] <- Read10X(data.dir = paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/Atamian/10x_files/", s))    
  Atamian_organoids [[s]] = CreateSeuratObject(counts = files_10x[[s]], project = paste0(s), min.cells = 3, min.features = 200)
}
#Rest of the code is adapted from official repository: https://github.com/quadratolab/Cerebellar-Organoid-2023/blob/main/integrateObject.v4.only.D2.manuscript.R
#QC - quite oversimplified in my opinion
for (obj in names(Atamian_organoids)) {
  Atamian_organoids[[obj]] = subset(x=Atamian_organoids[[obj]], cells=rownames(Atamian_organoids[[obj]]@meta.data[Atamian_organoids[[obj]]@meta.data$nFeature_RNA>800,]))
}
# Get the percentage mitochondrial genes
Atamian_organoids[[1]][["percent.mt"]] <- PercentageFeatureSet(Atamian_organoids[[1]], pattern = "^MT-")
Atamian_organoids[[2]][["percent.mt"]] <- PercentageFeatureSet(Atamian_organoids[[2]], pattern = "^MT-")
Atamian_organoids[[3]][["percent.mt"]] <- PercentageFeatureSet(Atamian_organoids[[3]], pattern = "^MT-")
# Set cutoff at 15% the right side of the normal distribution
Atamian_organoids[[1]]<- subset(Atamian_organoids[[1]], subset = percent.mt < 15)
Atamian_organoids[[2]]<- subset(Atamian_organoids[[2]], subset = percent.mt < 15)
Atamian_organoids[[3]]<- subset(Atamian_organoids[[3]], subset = percent.mt < 15)

#SCT Normalization
Atamian_organoids <- lapply(X = Atamian_organoids, FUN = SCTransform, ncells=5000, conserve.memory=T)
s.genes <- cc.genes$s.genes # automatically in Seurat package.
g2m.genes <- cc.genes$g2m.genes
#Calculate the CC phase and then regress them out
Atamian_organoids[[1]]<- CellCycleScoring(Atamian_organoids[[1]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Atamian_organoids[[2]]<- CellCycleScoring(Atamian_organoids[[2]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Atamian_organoids[[3]]<- CellCycleScoring(Atamian_organoids[[3]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Atamian_organoids[[1]][["CC.Difference"]] <- Atamian_organoids[[1]]$S.Score - Atamian_organoids[[1]]$G2M.Score
Atamian_organoids[[1]]$old.ident <- NULL
Atamian_organoids[[2]][["CC.Difference"]] <- Atamian_organoids[[2]]$S.Score - Atamian_organoids[[2]]$G2M.Score
Atamian_organoids[[2]]$old.ident <- NULL
Atamian_organoids[[3]][["CC.Difference"]] <- Atamian_organoids[[3]]$S.Score - Atamian_organoids[[3]]$G2M.Score
Atamian_organoids[[3]]$old.ident <- NULL
Atamian_organoids <- lapply(X = Atamian_organoids, FUN = SCTransform, ncells=5000, conserve.memory=T)
#Integration
all_features <- lapply(Atamian_organoids, row.names) %>% Reduce(intersect, .) 
features <- SelectIntegrationFeatures(object.list = Atamian_organoids, nfeatures = 3000)
Atamian_organoids <- PrepSCTIntegration(object.list = Atamian_organoids, anchor.features = features)
object.anchors <- FindIntegrationAnchors(object.list = Atamian_organoids, normalization.method = "SCT", anchor.features = features)
Atamian_integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT", features.to.integrate = all_features)
#DimensionalityReduction
Atamian_integrated <- RunPCA(Atamian_integrated, verbose = T)
Atamian_integrated <- RunUMAP(Atamian_integrated, reduction = "pca", dims = 1:30)
#Attach Metadata
Meta = read.delim("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/Atamian/Metadata/GSE247974_metadata.txt")
Clusters = as.data.frame(Meta$final.clusters, row_names = 1)
rownames(Clusters) = Meta$UMI
Atamian_integrated = AddMetaData(Atamian_integrated, Clusters, "final_clusters")
DimPlot(Atamian_integrated, group.by= "final_clusters", label = T)
Idents(Atamian_integrated) = "final_clusters"
# Subset to exclude cells with NA in the 'final_clusters'
Celltypes = unique(Atamian_integrated$final_clusters)
Celltypes = as.character(na.omit(Celltypes))
Atamian_integrated <- subset(Atamian_integrated, idents = Celltypes)

saveRDS(Atamian_integrated, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/Atamian/Seurat_files/Atamian_integrated.rds")
