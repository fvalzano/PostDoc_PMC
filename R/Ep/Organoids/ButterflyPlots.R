library(Seurat)
library(harmony)
library(SCpubr)
library(ggplot2)
library(tidyverse)
scrna_ha_pos = readRDS("/hpc/pmc_kool/mroosen/gebtoanalysis/RDS files/d30_d45_integrated_ha_only.rds")
scrna_ha_pos_list = SplitObject(scrna_ha_pos, split.by= "sample")
for (i in names(scrna_ha_pos_list)) {
  DefaultAssay(scrna_ha_pos_list[[i]]) = "RNA"
  scrna_ha_pos_list[[i]] = SCTransform(scrna_ha_pos_list[[i]], vars.to.regress = "pct_counts_mt", vst.flavor = "v2", assay = "RNA")
}
hvg= SelectIntegrationFeatures(scrna_ha_pos_list, nfeatures = 3000)
scrna_ha_pos = merge(x = scrna_ha_pos_list[[1]], y= scrna_ha_pos_list[-1], merge.data = TRUE, project = "Ep") 
scrna_ha_pos = RunPCA (scrna_ha_pos, verbose = FALSE, assay = "SCT", npcs= 50, features = hvg)
scrna_ha_pos = RunHarmony(scrna_ha_pos, group.by.vars = "sample", reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
scrna_ha_pos = RunUMAP(scrna_ha_pos, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
scrna_ha_pos = FindNeighbors(object = scrna_ha_pos, reduction = "harmony", dims = 1:30)
i = seq(0.1, 2, by = 0.1)
scrna_ha_pos = FindClusters(scrna_ha_pos, resolution = i)
#Load yap and export metadata 'own.mapping'
scrna_yap = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/HATag_Positive/YAP_tumor_20240322.rds")
own.mapping.yap = as.data.frame(scrna_yap$own_mapping)
colnames(own.mapping.yap) = "own_mapping"
#Load zfta and export metadata 'own.mapping'
scrna_zfta = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/HATag_Positive/zfta_tumor.rds")
own.mapping.zfta = as.data.frame(scrna_zfta$own_mapping)
colnames(own.mapping.zfta) = "own_mapping"
#Merge the two metadata
own.mapping.merge = rbind(own.mapping.yap, own.mapping.zfta)
#Add merge emtadata from yap and zfta onto the integrated yap-zfta object
scrna_ha_pos = AddMetaData(scrna_ha_pos, own.mapping.merge, col.name = "own_mapping")
#Split the integrated zfta-yap object
Idents(scrna_ha_pos) = "subtype"
scrna_ha_pos_yap = subset(scrna_ha_pos, idents = "YAP")
scrna_ha_pos_zfta = subset(scrna_ha_pos, idents = "ZFTA")
#Load the healthy organoids object to obtain healthy cluster-specific markers
scrna = readRDS("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/rds files/LX_healthy_final.rds")
Idents(scrna) = "own.mapping"
scrna = RenameIdents(scrna, 
                    c("Early RG" = "Early RG",
                      "Neuronal" = "Neuronal",
                      "Late RG" = "Late RG",
                      "RG/Astroglia" = "EPN/ECM"))
scrna$own.mapping = scrna@active.ident
DEG = FindAllMarkers(scrna, min.diff.pct=0.15, min.pct=0.1)
DEG = DEG[DEG$p_val_adj <= 0.05,]
DEG = DEG[order(DEG$avg_log2FC, decreasing = T),]
gene_list = list()
for(i in DEG$cluster){
    gene_list[[i]] = head(DEG[DEG$cluster == i,]$gene, n = 30)
}
#Avoid spaces in list names
names(gene_list) = c("Early_RG", "EPN/ECM", "Neuronal", "Late_RG")
#Export the plot for ggplotting
p = do_CellularStatesPlot(sample = scrna_ha_pos, 
                      input_gene_list= gene_list,
                      group.by ="subtype",
                      x1 = "Early_RG",
                      x2 = "Late_RG",
                      y1 = "EPN/ECM",
                      y2 = "Neuronal")
bplot = p$data
bplot$barcodes = rownames(bplot)
own.mapping.merge$barcodes = rownames(own.mapping.merge)
bplot_merge = merge(bplot, own.mapping.merge, by = "barcodes")
#Export the plots
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Ependymoma_organoids/Plots/ButterflyPlot.pdf", width = 7.5, height = 5)
ggplot(bplot_merge%>% arrange(group.by), aes(x = set_x, y = set_y))+
    geom_point(aes(colour = group.by))+
    theme_bw()+
    annotate("text", x = -Inf, y = Inf, label = "Early RG",
           hjust = -0.1, vjust = 1.1, size = 7.5, color = "black")+
    annotate("text", x = Inf, y = Inf, label = "Late RG",
           hjust = 1.1, vjust = 1.1, size = 7.5, color = "black")+
    annotate("text", x = -Inf, y = -Inf, label = "EPN/ECM",
           hjust = -0.1, vjust = -1.1, size = 7.5, color = "black")+
    annotate("text", x = Inf, y = -Inf, label = "Neuronal",
           hjust = 1.1, vjust = -1.1, size = 7.5, color = "black")+
    xlab("")+ylab("")+
    scale_colour_manual(values = c("YAP" = "#F6BD60", "ZFTA" = "#84A59D"))
ggplot(bplot_merge, aes(x = set_x, y = set_y))+
    geom_point(aes(colour = own_mapping))+
    theme_bw()+
    annotate("text", x = -Inf, y = Inf, label = "Early RG",
           hjust = -0.1, vjust = 1.1, size = 7.5, color = "black")+
    annotate("text", x = Inf, y = Inf, label = "Late RG",
           hjust = 1.1, vjust = 1.1, size = 7.5, color = "black")+
    annotate("text", x = -Inf, y = -Inf, label = "EPN/ECM",
           hjust = -0.1, vjust = -1.1, size = 7.5, color = "black")+
    annotate("text", x = Inf, y = -Inf, label = "Neuronal",
           hjust = 1.1, vjust = -1.1, size = 7.5, color = "black")+
    xlab("")+ylab("")+
    scale_colour_manual(values = c("EPN_ECM" = "#faac6f", 
                                   "Neuronal" = "#d9a27a",
                                   "Late_RG" = "#f0e8f2",
                                   "Early_RG" = "#6b4783",
                                   "mixed" = "#aaaaaa"))
dev.off()
