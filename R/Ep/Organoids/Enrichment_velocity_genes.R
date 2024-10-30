library(enrichR)
library(DOSE)
library(clusterProfiler)
library(enrichplot)

#Load velocity genes from timepoints single cell RNA sequencing experiment 
genes = read.csv("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/Timepoints/top_genes_ordered.csv")
genes = as.data.frame(t(genes))
genes$Gene_symbol = rownames(genes)
#Enrichment analysis
dbs = "GO_Biological_Process_2023"
Enrichment = enrichr(genes$Gene_symbol, dbs)
Enrichment_df = as.data.frame(Enrichment$GO_Biological_Process_2023)
Enrichment_df = Enrichment_df[order(Enrichment_df$Adjusted.P.value, decreasing = F),]
Enrichment_df = head(Enrichment_df, n = 20)
Genes = Enrichment_df$Genes
Genes = strsplit(Genes, ";")
Enrichment_df$Nr.Genes = c(lengths(Genes[1]),
                           lengths(Genes[2]),
                           lengths(Genes[3]),
                           lengths(Genes[4]),
                           lengths(Genes[5]),
                           lengths(Genes[6]),
                           lengths(Genes[7]),
                           lengths(Genes[8]),
                           lengths(Genes[9]),
                           lengths(Genes[10]),
                           lengths(Genes[11]),
                           lengths(Genes[12]),
                           lengths(Genes[13]),
                           lengths(Genes[14]),
                           lengths(Genes[15]),
                           lengths(Genes[16]),
                           lengths(Genes[17]),
                           lengths(Genes[18]),
                           lengths(Genes[19]),
                           lengths(Genes[20])
                           )
#Very basic bubble plot - maybe improve?
pdf("/hpc/pmc_kool/fvalzano/Jupyter/scvelo/Ependymoma/Timepoints/figures/Enrichment_results.pdf", width = 7.5, height = 5)
ggplot(Enrichment_df, aes(y = Enrichment_df$Term, x = -log10(Enrichment_df$Adjusted.P.value), size = Enrichment_df$Nr.Genes)) +
    geom_point()+
    xlab("Significance (-log10(Adjusted P values))") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 10))+
    guides(size=guide_legend(title="Nr.Genes"))
dev.off()
