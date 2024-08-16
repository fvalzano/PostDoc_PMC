library(Seurat)

scrna_mimic = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
Idents(scrna_mimic) = "SCT_snn_res.0.4"
scrna_mimic_immune = subset(scrna_mimic, idents = "10")
scrna_mimic_immune = RunPCA(scrna_mimic_immune, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_mimic_immune <- IntegrateLayers(object = scrna_mimic_immune, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT", normalization.method = "SCT",
verbose = FALSE)
scrna_mimic_immune = RunUMAP(scrna_mimic_immune, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
DimPlot(scrna_mimic_immune, pt.size = 1.25, group.by = "Cancer.type", reduction = "umap_harmony")


