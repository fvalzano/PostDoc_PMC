library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)

scrna_visium <- Load10X_Spatial(data.dir = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Data_examples/Visium_HD/Mouse_brain/10x/", bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(scrna_visium)
DefaultAssay(scrna_visium) <- "Spatial.008um"

options(vsc.dev.args = list(width=500, height=500, pointsize=6, res=100))
vln.plot <- VlnPlot(scrna_visium, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
#count.plot <- 
SpatialFeaturePlot(scrna_visium, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot

# Normalization both 8um and 16um bins - SCT is also available
DefaultAssay(scrna_visium) <- "Spatial.008um"
scrna_visium <- NormalizeData(scrna_visium)
DefaultAssay(scrna_visium) <- "Spatial.016um"
scrna_visium <- NormalizeData(scrna_visium)

# Test plotting
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Data_examples/Plots/Test_plots.pdf", width = 500, height = 500)
DefaultAssay(scrna_visium) <- "Spatial.016um"
SpatialFeaturePlot(scrna_visium, features = "Hpca") + ggtitle("Hpca expression (16um)")
# switch back to 8um
DefaultAssay(scrna_visium) <- "Spatial.008um"
SpatialFeaturePlot(scrna_visium, features = "Hpca") + ggtitle("Hpca expression (8um)")
dev.off()

# Unsupervised clustering
DefaultAssay(scrna_visium) <- "Spatial.008um"
scrna_visium <- FindVariableFeatures(scrna_visium)
scrna_visium <- ScaleData(scrna_visium)
# Sketch method with 50,0000 cells - This creates a new 'sketch' assay
scrna_visium <- SketchData(
  scrna_visium = scrna_visium,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
# switch analysis to sketched cells
DefaultAssay(scrna_visium) <- "sketch"

# perform clustering workflow
scrna_visium <- FindVariableFeatures(scrna_visium)
scrna_visium <- ScaleData(scrna_visium)
scrna_visium <- RunPCA(scrna_visium, assay = "sketch", reduction.name = "pca.sketch")
scrna_visium <- FindNeighbors(scrna_visium, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
scrna_visium <- FindClusters(scrna_visium, cluster.name = "seurat_cluster.sketched", resolution = 3)
scrna_visium <- RunUMAP(scrna_visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
DefaultAssay(scrna_visium) <- "sketch"
Idents(scrna_visium) <- "seurat_cluster.sketched"
DimPlot(scrna_visium, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") 
scrna_visium <- ProjectData(
  object = scrna_visium,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# switch to full dataset
DefaultAssay(scrna_visium) <- "Spatial.008um"
Idents(scrna_visium) <- "seurat_cluster.projected"
DimPlot(scrna_visium, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") 


SpatialDimPlot(scrna_visium, label = T, repel = T, label.size = 4)
Idents(scrna_visium) <- "seurat_cluster.projected"
cells <- CellsByIdentities(scrna_visium, idents = c(0, 4, 32, 34, 35))
p <- SpatialDimPlot(scrna_visium,
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p


#Identification of spatially-defined tissue domains with Banksy requires R 4.4 -> Install when we will need this