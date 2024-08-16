library(Seurat)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(viridis)
library(dittoSeq)
library(dplyr)
library(str2str)
library(presto)
scrna = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Integration/scrna_harmony.rds")
#Overview public datasets
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_FullDataset.pdf", width = 7.5, height = 5)
DimPlot(scrna,
        reduction = "umap",
        group.by = "Dataset",
        raster = F)+
        ggtitle("Unintegrated dataset (No batch-effect correction performed)")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
DimPlot(scrna,
        reduction = "umap_harmony",
        group.by = "Dataset",
        raster = F)+
        ggtitle("Integrated dataset (Harmony batch-effect correction)")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
dev.off()

#Overview public datasets with information regarding Entity and tumor subtypes
Idents(scrna) = "Entity"
scrna = subset(scrna, idents = c("Ependymoma", "Medulloblastoma"))
## Generate color list from Rcolorbrewer
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
## Split the color list in two sublist with distinct colours - Like this we avoid that the two plots will have same colors
col_list = list()
col_list[["Ependymoma"]] = unique(col_vector[1:length(col_vector)/2])
col_list[["Medulloblastoma"]] = unique(col_vector[length(col_vector)/2:length(col_vector)])
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_splitby_Entity.pdf", width = 10, height = 5)
DimPlot(scrna,
        reduction = "umap_harmony",
        split.by = "Entity",
        group.by = "Entity",
        raster = F, 
        cols = sample(col_vector, length(unique(scrna$Entity))))+
        ggtitle("")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_splitby_Entity_groupby_Subgroup.pdf", width = 7.5, height = 5)
for (i in unique(scrna$Entity)){
    scrna_sub = list()
    scrna_sub[[i]] = subset(scrna, idents = i)
    p = DimPlot(scrna_sub[[i]],
        reduction = "umap_harmony",
        split.by = "Entity",
        group.by = "Subgroup",
        raster = F, 
        cols = col_list[[i]])+
        ggtitle("")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
    print(p)
    rm(scrna_sub)
    print(paste0("Done ", i))
}
dev.off()

#Fulldataset with clustering resolution 0.4
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_FullDataset_Lou_Clustering.pdf", width = 7.5, height = 5)
DimPlot(scrna,
        reduction = "umap_harmony",
        group.by = "SCT_snn_res.0.4",
        raster = F,
        cols = col_vector[1:length(unique(scrna$SCT_snn_res.0.4))])+
        ggtitle("Louvain Clustering (Resolution: 0.4)")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
dev.off()
#Fulldataset with expression of the immune cell marker PTPRC(CD45)
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/FeaturePlot_PTPRC.pdf", width = 7.5, height = 5)
FeaturePlot(scrna,
        reduction = "umap_harmony",
        features = "PTPRC",
        raster = F,
        order = T) +
        ggtitle("PTPRC (aka CD45)")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
dev.off()

rm(scrna)
## Overview of the Immune data subset
scrna_immune = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Immune_Splitting/scrna_immune.rds")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_Immune_Lou_Clustering.pdf", width = 7.5, height = 5)
DimPlot(scrna_immune,
        reduction = "umap_harmony",
        group.by = "SCT_snn_res.0.4",
        raster = F)+
        ggtitle("Louvain Clustering (Resolution: 0.4)")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
dev.off()
#Create a pseudobulk of the immune cell subcluster and perform correlation of the pseudobulk with lymphoid signature (CD3D) to identify lymphoid (High) and Myeloid (Low) compartments
pseudo = AggregateExpression(scrna_immune, 
                                      assays = "SCT", 
                                      return.seurat = T,
                                      group.by = "SCT_snn_res.0.4")
pseudo$SCT_snn_res.0.4 = paste(scrna_immune$SCT_snn_res.0.4)
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/Immune_correlation.pdf", width = 5, height = 5)
FeatureScatter(pseudo, 
               feature1="CD3D", 
               feature2="ITGAX",
               pt.size = 5)+
        ggtitle("PseudoBulk correlation")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
dev.off()

# Overview of Myeloid and Lymphoid compartment (UMAP + Dotplot)
scrna_immune_myeloid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Annotation/scrna_immune_myeloid.rds")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_Myeloid_Anno_v2.pdf", width = 7.5, height = 5)
DimPlot(scrna_immune_myeloid,
        reduction = "umap_harmony",
        group.by = "annotation_fv_v2",
        raster = F)+
        ggtitle("Myeloid cell compartment: annotation version v2")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
dev.off()
Idents(scrna_immune_myeloid) = "annotation_fv_v2"
#Retrieve the top n (in this case 3 for ease of view)
myeloid_markers = FindAllMarkers(scrna_immune_myeloid, 
                                      only.pos = T, 
                                      min.pct = 0.15,
                                      min.diff.pct = 0.15)
top_markers = list()
for (i in unique(myeloid_markers$cluster)) {
   markers = myeloid_markers[myeloid_markers$cluster == i,]
   markers = markers[markers$p_val_adj <= 0.05,]
   markers = markers[order(markers$avg_log2FC, decreasing = T),]
   top_markers[[i]] = head(markers, n = 3)   
}
top_markers = Join(data.list = top_markers, by = "gene")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/Dotplot_Myeloid_top3.pdf", width = 15, height = 5)
DotPlot(scrna_immune_myeloid,
        features = top_markers$gene,
        col.min = 0)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_line(linewidth = 0.1))+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
dev.off()

scrna_immune_lymphoid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Annotation/scrna_immune_lymphoid.rds")
scrna_immune_lymphoid = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Annotation/scrna_immune_lymphoid.rds")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/UMAP_Lymphoid_Anno_v2.pdf", width = 7.5, height = 5)
DimPlot(scrna_immune_lymphoid,
        reduction = "umap_harmony",
        group.by = "annotation_fv_v2",
        raster = F)+
        ggtitle("Lymphoid cell compartment: annotation version v2")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
dev.off()
Idents(scrna_immune_lymphoid) = "annotation_fv_v2"
#Retrieve the top n (in this case 3 for ease of view)
lymphoid_markers = FindAllMarkers(scrna_immune_lymphoid, 
                                      only.pos = T, 
                                      min.pct = 0.1,
                                      min.diff.pct = 0.1)
top_markers = list()
for (i in unique(lymphoid_markers$cluster)) {
   markers = lymphoid_markers[lymphoid_markers$cluster == i,]
   markers = markers[markers$p_val_adj <= 0.05,]
   markers = markers[order(markers$avg_log2FC, decreasing = T),]
   top_markers[[i]] = head(markers, n = 3)   
}
top_markers = Join(data.list = top_markers, by = "gene")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/Dotplot_Lymphoid_top3.pdf", width = 15, height = 5)
DotPlot(scrna_immune_lymphoid,
        features = top_markers$gene,
        col.min = 0)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_line(linewidth = 0.1))+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
dev.off()

##MIMIC section
scrna_mimic = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/MIMIC_UMAP.pdf", width = 7.5, height = 5)
DimPlot(scrna_mimic, reduction = "umap", group.by="Cancer.type") +        
        ggtitle("Unintegrated Dataset (No batch effect applied)")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
DimPlot(scrna_mimic, reduction = "umap_harmony", group.by="Cancer.type")+
        ggtitle("Integrated Dataset (Harmony batch effect removal applied)")+
        theme(title = element_text(size = 10), legend.text = element_text(size = 10))
dev.off()
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/MIMIC_UMAP_Lou0.4.pdf", width = 7.5, height = 5)
DimPlot(scrna_mimic, reduction = "umap_harmony", group.by= "SCT_snn_res.0.4", label = T)+
        ggtitle("Louvain Clustering (Resolution: 0.4)")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))
scrna_mimic$Major_classes_FV = ifelse(scrna_mimic$SCT_snn_res.0.4 == "10", "Immune",
                                        ifelse(scrna_mimic$SCT_snn_res.0.4 == "21", "Astroglia",
                                                ifelse(scrna_mimic$SCT_snn_res.0.4 == "12", "Fb1",
                                                        ifelse(scrna_mimic$SCT_snn_res.0.4 == "23", "Fb2",
                                                                ifelse(scrna_mimic$SCT_snn_res.0.4 == "30", "Pericytes",
                                                                        ifelse(scrna_mimic$SCT_snn_res.0.4 == "24", "EC", "Neuronal"))))))
DimPlot(scrna_mimic, group.by= "Major_classes_FV", reduction = "umap_harmony")+
        ggtitle("Major cell lineages")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))

dev.off()

pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/MIMIC_UMAP_Immune_Subset.pdf", width = 7.5, height = 5)
scrna_mimic = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
Idents(scrna_mimic) = "SCT_snn_res.0.4"
scrna_mimic_immune = subset(scrna_mimic, idents = "10")
scrna_mimic_immune = RunPCA(scrna_mimic_immune, verbose = FALSE, assay = "SCT", npcs= 50)
scrna_mimic_immune <- IntegrateLayers(object = scrna_mimic_immune, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT", normalization.method = "SCT",
verbose = FALSE)
scrna_mimic_immune = RunUMAP(scrna_mimic_immune, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
DimPlot(scrna_mimic_immune, group.by = "Cancer.type", reduction = "umap_harmony")+
        ggtitle("Immune Cells subset")+
        theme(title = element_text(size = 10), 
              legend.text = element_text(size = 10))

dev.off()
pdf("/hpc/pmc_kool/fvalzano/Presentation_plots/20240822_TME_update/MIMIC_GeneExpression_Immune_Subset.pdf", width = 7.5, height = 5)
FeaturePlot(scrna_mimic_immune, features = "ITGAX", order = T, reduction = "umap_harmony")+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
FeaturePlot(scrna_mimic_immune, features = "CD3D", order = T, reduction = "umap_harmony")+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
FeaturePlot(scrna_mimic_immune, features = "CD19", order = T, reduction = "umap_harmony")+
        scale_color_gradientn(colors = plasma(n = 10, direction = -1))
dev.off()