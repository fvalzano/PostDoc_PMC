library(Seurat)
library(readxl)
library(ggplot2)
library(circlize)
library(readr)
#Load the healthy organoid dataset
scrna = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")

#Load the NMF metaprogram from Filbin et al
NMF_meta = read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Annotations/NMFmetaprorams_from_Filbin/media-5-2.xlsx")
colnames(NMF_meta) = NMF_meta[1,]
NMF_meta[1,] = NA
NMF_meta = na.omit(NMF_meta)
for(i in colnames(NMF_meta)){
    scrna = AddModuleScore(scrna, features = list(NMF_meta[[i]]), name = paste0(i, "_score"))
}
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/Overlay_Filbin_unbiased_clustering.pdf", width = 10, height = 7.5)
DimPlot(scrna,
    group.by = "RNA_snn_res.0.5", 
    label = T) +
DotPlot(scrna, features = c("Neuroepithelial-like_score1",
                            "MES/Hypoxia_score1",
                            "Ependymal-like_score1",
                            "Radial-glia-like_score1",
                            "Cycling_score1",
                            "Embryonic-neuronal-like_score1",
                            "Neuronal-like_score1",
                            "Embryonic-like_score1"), 
               col.min = 0, 
               group.by = "RNA_snn_res.0.5")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/Overlay_Filbin_ownmapping.pdf", width = 10, height = 7.5)
DimPlot(scrna,
    group.by = "own.mapping", 
    label = T) +
DotPlot(scrna, features = c("Neuroepithelial-like_score1",
                            "MES/Hypoxia_score1",
                            "Ependymal-like_score1",
                            "Radial-glia-like_score1",
                            "Cycling_score1",
                            "Embryonic-neuronal-like_score1",
                            "Neuronal-like_score1",
                            "Embryonic-like_score1"), 
               col.min = 0, 
               group.by = "own.mapping")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

p = DotPlot(scrna, features = c("Neuroepithelial-like_score1",
                            "MES/Hypoxia_score1",
                            "Ependymal-like_score1",
                            "Radial-glia-like_score1",
                            "Cycling_score1",
                            "Embryonic-neuronal-like_score1",
                            "Neuronal-like_score1",
                            "Embryonic-like_score1"), 
               col.min = 0, 
               group.by = "own.mapping")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


chord = data.frame("Filbin.annotation" = c("Neuroepithelial-like",
                            "MES/Hypoxia",
                            "Ependymal-like",
                            "Radial-glia-like",
                            "Cycling",
                            "Cycling",
                            "Embryonic-neuronal-like",
                            "Neuronal-like",
                            "Embryonic-like" 
                                           ),
                    "Own.mapping" = c("RG/Astroglia",
                                      "RG/Astroglia",
                                      "Neuronal",
                                      "Neuronal",
                                      "Late RG",
                                      "Early RG",
                                      "Neuronal",
                                      "Neuronal",
                                      "Neuronal"))

pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/Chord_Ownmapp_Filbin.pdf", width = 5, height = 5)
chordDiagram(chord, 
             annotationTrack = c("name", "grid"),
             grid.col = c("RG/Astroglia" = "#F8766D",
                          "Neuronal" = "#7CAE00",
                          "Late RG" = "#00BFC4",
                          "Early RG" = "#C77CFF"))
dev.off()

write_rds(scrna, "/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
