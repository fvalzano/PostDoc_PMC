.libPaths("/hpc/local/Rocky8/pmc_kool/R-4.1.2/lib64/R/x86_64-pc-linux-gnu-library/4.1")
library(readr)
library(Seurat)
library(SeuratDisk)


whole_zfta = read_rds("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/whole_zfta-1.rds")
whole_zfta@meta.data[,5:25] =NULL
whole_zfta$own_mapping = as.character(whole_zfta$own_mapping)
whole_zfta@assays$predictions.linnarson = NULL
whole_zfta@assays$RNA$data = NULL
whole_zfta@assays$RNA$scale.data = NULL
#Seurat v5 altered the conversion performed by SeuratDisk -> Convert RNA assay to v3-like format
whole_zfta[["RNA3"]] <- as(whole_zfta[["RNA"]], Class = "Assay")
DefaultAssay(whole_zfta) <- "RNA3"
whole_zfta[["RNA"]] <- NULL
whole_zfta <- RenameAssays(whole_zfta, RNA3 = 'RNA')
SaveH5Seurat(whole_zfta, filename="/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/whole_zfta.h5Seurat", overwrite=T)
Convert("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/whole_zfta.h5Seurat", dest = "h5ad", overwrite=T)

zfta_wam = read_rds("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/zfta_wam.rds")
zfta_wam@assays$predictions.linnarson = NULL
zfta_wam@assays$RNA$data = NULL
zfta_wam@assays$RNA$scale.data = NULL
zfta_wam[["RNA3"]] <- as(zfta_wam[["RNA"]], Class = "Assay")
DefaultAssay(zfta_wam) <- "RNA3"
zfta_wam[["RNA"]] <- NULL
zfta_wam <- RenameAssays(zfta_wam, RNA3 = 'RNA')
SaveH5Seurat(zfta_wam, filename="/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/zfta_wam.h5Seurat", overwrite=T)
Convert("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/zfta_wam.h5Seurat", dest = "h5ad", overwrite=T)

scrna = read_rds("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
scrna@assays$predictions.linnarson = NULL
scrna@assays$RNA$data = NULL
scrna@assays$RNA$scale.data = NULL
scrna[["RNA3"]] <- as(scrna[["RNA"]], Class = "Assay")
DefaultAssay(scrna) <- "RNA3"
scrna[["RNA"]] <- NULL
scrna <- RenameAssays(scrna, RNA3 = 'RNA')
Idents(scrna) = "sample"
scrna = RenameIdents(scrna, c("LX146_d11" = "d11",
                              "LX147_d11" = "d11",
                              "LX117_d30" = "d30",
                              "LX362_d30" = "d30",
                              "LX370_d44" = "d44",
                              "LX086_d45" = "d44",
                              "LX138_d80" = "d80",
                              "LX141_d80" = "d80"))
scrna$timepoint = scrna@active.ident 
DimPlot(scrna, group.by = "timepoint")
Idents(scrna) = "timepoint"
scrna_list = list()
for (i in unique(Idents(scrna))) {
     scrna_list[[i]] = subset(scrna, idents = i)
     SaveH5Seurat(scrna_list[[i]], filename=paste0("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/", "scrna", i, ".h5Seurat"), overwrite=T)
     Convert(paste0("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/", "scrna", i, ".h5Seurat"), dest = "h5ad", overwrite=T)
}

SaveH5Seurat(scrna, filename="/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/scrna_timepoints.h5Seurat", overwrite=T)
Convert("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/scrna_timepoints.h5Seurat", dest = "h5ad", overwrite=T)
