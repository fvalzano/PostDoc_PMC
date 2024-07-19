# packages required for Visium HD

install.packages("hdf5r")
install.packages("arrow")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

object <- Load10X_Spatial(data.dir = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Data_examples/Visium_HD/Mouse_brain/10x/", bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot