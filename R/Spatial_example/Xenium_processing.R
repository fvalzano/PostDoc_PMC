#DISCLAIMER_FV: perform Xenium analysis on Jobhopper, as VSCode instances have troubles with dependencies of 'sf' R package,
#necessary for the cropping analysis. One of these dependencies is as 'units' <--- Fix this someday


library(Seurat)
library(future)
#plan("multisession", workers = 10)
library(ggplot2)

data_dir <- "/hpc/pmc_kool/fvalzano/Rstudio_Test1/Data_examples/Xenium"
# Load the Xenium data
xenium.obj <- LoadXenium(data_dir, fov = "fov")
# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
#Visualize data
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)
ImageFeaturePlot(xenium.obj, features = c("Cux2", "Rorb"), max.cutoff = c(25,
    35), size = 0.75, cols = c("white", "red"))

#PROBLEM: Missing dependencies
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
    coord.fixed = FALSE, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1"), nmols = 10000)