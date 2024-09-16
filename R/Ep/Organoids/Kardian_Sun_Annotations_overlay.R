library(readxl)
library(reshape2)
library(stringr)
library(Seurat)
library(ggplot2)
Hua_anno = read_xlsx("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Annotations/Annotation_from_Kardian_Sun/TableS2.human.xlsx")
#Sort phenotype needed for the annotation
Hua_anno = Hua_anno[Hua_anno$'Short name' %in% c("RGC", "nIPC", "Astro", "Neuron", "OPC"),]
#Split the gene column in single genes for phenotype lists
Markers = str_split(Hua_anno$Markers, ",")
#Name the lists according to the phenotypes
names(Markers) = Hua_anno$'Short name'
###Subset RGC phenotype###
RGC_like_genes = Markers[names(Markers) %in% "RGC"]
#Convert the lists into character
RGC_like_geneslist = character()
RGC_like_geneslist = RGC_like_genes[[1]]
for(i in names(RGC_like_genes)){
    RGC_like_genes[[1]] = NULL
    RGC_like_geneslist = c(RGC_like_geneslist, RGC_like_genes[[i]])
}
#Get rid of duplicated
RGC_like_geneslist = RGC_like_geneslist[duplicated(RGC_like_geneslist)]

###Subset nIPC phenotype###
nIPC_like_genes = Markers[names(Markers) %in% "nIPC"]
#Convert the lists into character
nIPC_like_geneslist = character()
nIPC_like_geneslist = nIPC_like_genes[[1]]
for(i in names(nIPC_like_genes)){
    nIPC_like_genes[[1]] = NULL
    nIPC_like_geneslist = c(nIPC_like_geneslist, nIPC_like_genes[[i]])
}
#Get rid of duplicated
nIPC_like_geneslist = nIPC_like_geneslist[duplicated(nIPC_like_geneslist)]

###Subset Astro phenotype###
Astro_like_genes = Markers[names(Markers) %in% "Astro"]
#Convert the lists into character
Astro_like_geneslist = character()
Astro_like_geneslist = Astro_like_genes[[1]]
for(i in names(Astro_like_genes)){
    Astro_like_genes[[1]] = NULL
    Astro_like_geneslist = c(Astro_like_geneslist, Astro_like_genes[[i]])
}
#Get rid of duplicated
Astro_like_geneslist = Astro_like_geneslist[duplicated(Astro_like_geneslist)]

###Subset Neuron phenotype###
Neuron_like_genes = Markers[names(Markers) %in% "Neuron"]
#Convert the lists into character
Neuron_like_geneslist = character()
Neuron_like_geneslist = Neuron_like_genes[[1]]
for(i in names(Neuron_like_genes)){
    Neuron_like_genes[[1]] = NULL
    Neuron_like_geneslist = c(Neuron_like_geneslist, Neuron_like_genes[[i]])
}
#Get rid of duplicated
Neuron_like_geneslist = Neuron_like_geneslist[duplicated(Neuron_like_geneslist)]

###Subset OPC phenotype###
OPC_like_genes = Markers[names(Markers) %in% "OPC"]
#Convert the lists into character
OPC_like_geneslist = character()
OPC_like_geneslist = OPC_like_genes[[1]]
for(i in names(OPC_like_genes)){
    OPC_like_genes[[1]] = NULL
    OPC_like_geneslist = c(OPC_like_geneslist, OPC_like_genes[[i]])
}
#Get rid of duplicated
OPC_like_geneslist = OPC_like_geneslist[duplicated(OPC_like_geneslist)]


#Load Miekanoids scrnaseq 
scrna = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
scrna = AddModuleScore(scrna, features = list(RGC_like_geneslist), name = "RGC_score")
scrna = AddModuleScore(scrna, features = list(nIPC_like_geneslist), name = "nIPC_score")
scrna = AddModuleScore(scrna, features = list(Astro_like_geneslist), name = "Astro_score")
scrna = AddModuleScore(scrna, features = list(Neuron_like_geneslist), name = "Neuron_score")
scrna = AddModuleScore(scrna, features = list(OPC_like_geneslist), name = "OPC_score")
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/Overlay_Kardian_Sun_ownmapping.pdf", width = 10, height = 7.5)
DimPlot(scrna,
    group.by = "own.mapping", 
    label = T) +
DotPlot(scrna, features = c("RGC_score1",
                            "nIPC_score1",
                            "Astro_score1",
                            "Neuron_score1",
                            "OPC_score1"), 
               col.min = 0, 
               group.by = "own.mapping")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

write_rds(scrna, "/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
