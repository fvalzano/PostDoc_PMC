library(Seurat)
library(ggplot2)
library(ggalluvial)
library(readxl)
library(stringr)
#Build dataframe with Miekanoids, Filbin and Mack annotations:
scrna = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
Idents(scrna) = "own.mapping"
scrna = RenameIdents(scrna, 
                    c("Early RG" = "Early RG",
                      "Neuronal" = "Neuronal",
                      "Late RG" = "Late RG",
                      "RG/Astroglia" = "EPN/ECM"))
scrna$own.mapping = scrna@active.ident
#Build dotplot based object to retrieve info about Filbin annotation pairs and expression
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
#The slot data contains the info we want
Miekanoids_Filbin = p$data
#Build dotplot based object to retrieve info about Kardian-Sun annotation pairs and expression
p = DotPlot(scrna, features = c("RGC_score1",
                            "nIPC_score1",
                            "Astro_score1",
                            "Neuron_score1",
                            "OPC_score1"), 
               col.min = 0, 
               group.by = "own.mapping")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
Miekanoids_KS= p$data
#Merge the information regarding annotation pairs and avg.scaled.expression 
Sankey1 = data.frame(Miekanoids_Filbin$features.plot, Miekanoids_Filbin$id, Miekanoids_Filbin$avg.exp.scaled)
colnames(Sankey1) = c("Annotation", "anno_mieke", "value")
#remove annotation pairs with 0 value
Sankey1$value = ifelse(Sankey1$value == 0, NA, Sankey1$value)
Sankey1 = na.omit(Sankey1)
Sankey1 = melt(Sankey1)
Sankey1$variable = NULL
Sankey2 = data.frame(Miekanoids_KS$features.plot, Miekanoids_KS$id, Miekanoids_KS$avg.exp.scaled)
colnames(Sankey2) = c("Annotation", "anno_mieke", "value")
#remove annotation pairs with 0 value
Sankey2$value = ifelse(Sankey2$value == 0, NA, Sankey2$value)
Sankey2 = na.omit(Sankey2)
Sankey2=melt(Sankey2)
Sankey2$variable=NULL
#merge the two Sankey df to create a unique df
Sankey = merge(Sankey1, Sankey2, by = "anno_mieke")
colnames(Sankey) = c("anno_mieke", "anno_filbin", "value_filbin", "anno_ks", "value_ks")
Sankey_melt = melt(Sankey)
Sankey_melt$anno_filbin = str_remove(Sankey_melt$anno_filbin, "_score1")
Sankey_melt$anno_ks = str_remove(Sankey_melt$anno_ks, "_score1")
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/Sankey_alluvial_merged_annotations.pdf", width = 10, height =10)
ggplot(Sankey_melt,
       aes(axis1 = anno_filbin, axis2 = anno_mieke , axis3 = anno_ks, y = value))+
       scale_x_discrete(limits = c("anno_filbin", "anno_mieke", "anno_ks"), labels = c("Filbin \nLab Annotation", "This \nStudy Annotation", "Mack \nLab Annotation"))+
       geom_alluvium(aes(fill=anno_mieke))+
       geom_stratum()+
       geom_text(stat="stratum", aes(label = after_stat(stratum)))+
       theme_minimal()+
       theme(panel.grid.major = element_blank(),
             axis.text.y = element_blank(),
             axis.title.y=element_blank(),
             axis.text.x=element_text(size=15, angle= 45, hjust = 1, vjust=1.25))
dev.off()
