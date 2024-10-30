library(enrichR)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
Plotting_directory = ""
#Load the list of significant DEG per comparison, we will focus for now on all the single comparisons against empty-cag as well as MDG vs MD plasmids
MYCN_vs_EMPTY_signif = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_vs_EMPTY.csv")
MYCN_DNTP53_vs_EMPTY_signif = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_vs_EMPTY.csv")
MYCN_DNTP53_GLI2_vs_EMPTY_signif = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_GLI2_vs_EMPTY.csv")
MYCN_DNTP53_GLI2_vs_MYCN.DNTP53_signif = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_GLI2_vs_MYCN.DNTP53.csv")
#Load all the databases from EnrichR
dbs_tot <- listEnrichrDbs()
#Focus on "GO_Biological_Process_2023", "GO_Molecular_Function_2023","KEGG_2021_Human" and "DisGeNET"
dbs = c("GO_Biological_Process_2023", 
        "GO_Molecular_Function_2023",
        "KEGG_2021_Human",
        "DisGeNET")
#Upregulated genes in MDG vs EMPTY
MDGvsEM = enrichr(MYCN_DNTP53_GLI2_vs_EMPTY_signif[MYCN_DNTP53_GLI2_vs_EMPTY_signif$log2FoldChange>1,]$X, dbs)
#Upregulated genes in MD vs EMPTY
MDvsEM = enrichr(MYCN_DNTP53_vs_EMPTY_signif[MYCN_DNTP53_vs_EMPTY_signif$log2FoldChange>1,]$X, dbs)
#Upregulated genes in M vs EMPTY
MvsEM = enrichr(MYCN_vs_EMPTY_signif[MYCN_vs_EMPTY_signif$log2FoldChange>1,]$X, dbs)
#Upregulated genes in MDG vs MD
MDG_vs_MD = enrichr(MYCN_DNTP53_GLI2_vs_MYCN.DNTP53_signif[MYCN_DNTP53_GLI2_vs_MYCN.DNTP53_signif$log2FoldChange>1,]$X, dbs)
#Upregulated genes in MD vs MDG
MD_vs_MDG = enrichr(MYCN_DNTP53_GLI2_vs_MYCN.DNTP53_signif[MYCN_DNTP53_GLI2_vs_MYCN.DNTP53_signif$log2FoldChange<1,]$X, dbs)

#The enrichment analysis performed on MDG_vs_MD with DisGeNET gave Childhood MB as significant result - maybe the triple combination is leading to a more
#MB like progression in the chrisganoids? 
MDG_vs_MD_DGN_top20 = head(MDG_vs_MD$DisGeNET, n = 10)
MDG_vs_MD_DGN_top20 = MDG_vs_MD_DGN_top20[order(MDG_vs_MD_DGN_top20$Adjusted.P.value, decreasing = F),]
#Plot the Enrichment results as a barplot
pdf(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/Plots/", Plotting_directory,"Enrichment_barplot_MDGvsMD_DisGeNET.pdf"), height= 7.5, width = 12.5)
ggplot(MDG_vs_MD_DGN_top20, aes(x = -log10(MDG_vs_MD_DGN_top20[,"Adjusted.P.value"]), y = MDG_vs_MD_DGN_top20[,"Term"], fill = MDG_vs_MD_DGN_top20$Adjusted.P.value))+
        geom_bar(stat="identity")+
        scale_y_discrete(limits = rev(MDG_vs_MD_DGN_top20$Term))+
        scale_fill_gradient(high="gray", low="#41598b")+
        geom_vline(xintercept = 1.3, linetype = "dashed")+
        theme_bw()+
        xlab("Significance")+
        ylab("DisGeNET_Term")+
        theme(axis.text.x=element_text(size = 10),
              axis.text.y=element_text(size = 15),
              axis.title = element_text(size = 20))
dev.off()

#Plot the enriched genes coming from DisGeNET and check how are they expressed in the single bulkRNA runs
#Filter the enrichment results for Childhood Medulloblastoma and get enriched genes
MB_enriched_genes = MDG_vs_MD_DGN_top20[MDG_vs_MD_DGN_top20$Term =="Childhood Medulloblastoma",]$Genes
#Split the chr vector in several elements - the separator is ;
MB_enriched_genes = strsplit(MB_enriched_genes, ";")
MB_enriched_genes = MB_enriched_genes[[1]]
#Load the vst normalized merged counts  
Bulk_RNA_vst = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/Bulk_RNA_merge_vst_Normalized.csv")
#Subset the general counts for the genes of interest(GOI)
Bulk_RNA_vst_subset = Bulk_RNA_vst[Bulk_RNA_vst$X %in% MB_enriched_genes,]
rownames(Bulk_RNA_vst_subset) = Bulk_RNA_vst_subset$X
Bulk_RNA_vst_subset$X = NULL
#Plot GOI as boxplot in the single conditions
for(i in rownames(Bulk_RNA_vst_subset)){
        Bulk_RNA_vst_subset_genes = Bulk_RNA_vst_subset[rownames(Bulk_RNA_vst_subset) %in% i,]
        Bulk_RNA_vst_subset_genes = reshape2::melt(Bulk_RNA_vst_subset_genes)
        Bulk_RNA_vst_subset_genes$grouping = c("MYCN", 
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY",
                                        "MYCN", 
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY",
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY")
        pdf(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/Plots/",Plotting_directory,i, "_Boxplot.pdf"), width = 5, height = 5)
        p = ggplot(Bulk_RNA_vst_subset_genes, aes(x = Bulk_RNA_vst_subset_genes$grouping, y = Bulk_RNA_vst_subset_genes$value))+
                geom_boxplot(aes(fill = Bulk_RNA_vst_subset_genes$grouping))+
                geom_jitter(width = 0, size = 2.5)+
                theme_bw()+
                ggtitle(paste0(i, " Epression"))+
                ylab("Normalized counts")+
                xlab("")+
                theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                        axis.text.y = element_text(size = 10),
                        title = element_text(size = 10),
                        legend.title=element_blank())
        print(p)
        dev.off()
}

#Plot GOI as connected scatterplot in the single conditions
Bulk_RNA_vst_subset_genes = Bulk_RNA_vst_subset[rownames(Bulk_RNA_vst_subset) %in% rownames(Bulk_RNA_vst_subset),]
Bulk_RNA_vst_subset_genes = reshape2::melt(Bulk_RNA_vst_subset_genes)
Bulk_RNA_vst_subset_genes$grouping = rep(c("MYCN", 
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY",
                                        "MYCN", 
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY",
                                        "MYCN_DNTP53", 
                                        "MYCN_DNTP53_GLI2", 
                                        "EMPTY"), each = 8)
Bulk_RNA_vst_subset_genes$genes = rep(c("NTN1", "GLI2", "BMP7", "CDK6", "PTCH2", "ERBB2", "EGFR", "PTCH1"))
average_counts <- Bulk_RNA_vst_subset_genes %>%
  group_by(grouping, genes) %>%
  summarize(average_value = median(value, na.rm = TRUE)) %>%
  ungroup()

pdf(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/Plots/",Plotting_directory,"Connected_scatter.pdf"), width = 5, height = 5)
ggplot(average_counts, aes(x = average_counts$grouping, y = average_counts$average_value))+
        geom_point(aes(size = 1, colour = factor(average_counts$genes)))+
        scale_size(c(1,5))+
        geom_line(aes(group = average_counts$genes, colour = factor(average_counts$genes)))+
        theme_bw() +
        xlab("")+
        ylab("Median Normalized Counts")+
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                        axis.text.y = element_text(size = 10),
                        title = element_text(size = 10))
dev.off()
#Focus on the MDG comparison - are genes related to specific pathways in developing CB reflected in this setting?
Bulk_RNA_vst = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/Bulk_RNA_merge_vst_Normalized.csv")
#List of genes from PMC10924233 figure 5
Developmental_genes = c("PAX6","ZIC1","ZIC2","PCNA","MK167","DLGAPS","CCND1","NEUROD1","NEUROD2","CNTN1","CNTN2","UNCSC","EOMES","CALB2","GRM1","TBR1","GLI1","GLI2","GLI3","PTCH1")
#Subset the general counts for the genes of interest(GOI)
Bulk_RNA_vst_subset = Bulk_RNA_vst[Bulk_RNA_vst$X %in% Developmental_genes,]
rownames(Bulk_RNA_vst_subset) = Bulk_RNA_vst_subset$X
Bulk_RNA_vst_subset$X = NULL
#Subset the sample of interest
Bulk_RNA_vst_subset = Bulk_RNA_vst_subset[,c("CB2402.03.mix.6.MYCN.DNTP53.GLI2","CB2402.03.mix.8.empty.cag.ig","CB2410.11.mix.6.MYCN.DNTP53.GLI2","CB2410.11.mix.8.empty.cag.ig","CB2412.13.mix.6.MYCN.DNTP53.GLI2","CB2412.13.mix.8.empty.cag.ig"),]
#Generate annotation dataframe for heatmap
anno_row=data.frame(c("Early and proliferating GCP",
                      "Post-mitotic and migrating GCP",
                      "Early and proliferating GCP",
                      "Sonic hedgehog signalling",
                      "Sonic hedgehog signalling",
                      "Early and proliferating GCP",
                      "Sonic hedgehog signalling",
                      "Early and proliferating GCP",
                      "Deep cerebellar nuclei neurons",
                      "Unipolar brush cells",
                      "Early and proliferating GCP",
                      "Post-mitotic and migrating GCP",
                      "Unipolar brush cells",
                      "Post-mitotic and migrating GCP",
                      "Unipolar brush cells",
                      "Post-mitotic and migrating GCP",
                      "Sonic hedgehog signalling"))
#Adapt row and colnames
rownames(anno_row) = rownames(Bulk_RNA_vst_subset)
colnames(anno_row) = "annotation"
pdf("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/Plots/Heatmap_Cerebellar_Genes.pdf", width = 7.5, height = 7.5)
pheatmap(Bulk_RNA_vst_subset, 
         scale= "row", 
         cluster_rows=T, 
         cluster_cols=T,
         annotation_row = anno_row,
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()
