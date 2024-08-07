library(DESeq2)
library(purrr)
library(apeglm)
library(ggplot2)
library(umap)
library(uwot)
#Load bulk RNAseq runs from specific requests folder
##----------------------IMPORTANT=Modify request folder variable and IDs----------------------
result_directory = "20240807_Julie"
##Set up list to contain the single bulk RNA seq runs, delete unnecessary columns and rename remaining ones (counts and gene name)
RNAseq_files = list.files(paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests/", result_directory))
RNAseq_runs = list()
for (file in RNAseq_files) {
   RNAseq_runs[[file]]=read.table(paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests/", result_directory, "/", file), )
   #Save only V2 and V11 as they store counts and gene_names respectively
   RNAseq_runs[[file]]= RNAseq_runs[[file]][, c("V2", "V11")]
   colnames(RNAseq_runs[[file]]) = c("Counts", "Gene_name")
}
head(RNAseq_runs$'PMLBM000AFK_PMCRZ292YIC_RNA-Seq.gene_id.exon.counts.txt')
##Merge the single runs in one
RNAseq_merged = RNAseq_runs[[1]]
for (i in 2:length(RNAseq_runs)) {
  RNAseq_merged <- cbind(RNAseq_merged, RNAseq_runs[[i]])
}

##Remove Gene_name columns duplicates
RNAseq_merged = as.data.frame(RNAseq_merged, row.names = RNAseq_merged$Gene_name)
RNAseq_merged = RNAseq_merged[, -grep("Gene_name", names(RNAseq_merged))]

##Rename dataframe with IDs 
colnames(RNAseq_merged) = c("541AAO", "235AAR", "568AAR_Primary", "693AAR", "846AAS", "872AAS", "244AAT", "423AAT", "568AAR_Relapse", "235AAU")
##Write analysis design and set up dds object
#Create a design vector with High or Low expression of Rloops based on order of the RNAseq_merged dataframe
design = as.data.frame(c("Low", "High", "Low", "High", "Low", "High", "High", "High", "Low", "High"))
colnames(design) = "R_loops_expression"
rownames(design) = colnames(RNAseq_merged)
design$IDs = rownames(design)
dds <- DESeqDataSetFromMatrix(countData = RNAseq_merged,
                              colData = design,
                              design = ~ R_loops_expression)

##QC the dds object
smallestGroupSize <- 3
retain <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[retain,]
##Differentially Expressed Gene Analysis with apeglm shrinkage of fold changes
DES_dds <- DESeq(dds)
res <- results(DES_dds)
resultsNames(DES_dds)
resLFC <- lfcShrink(DES_dds, coef="R_loops_expression_Low_vs_High", type="ashr")
DEG = as.data.frame(resLFC)
DEG_significant = DEG[DEG$padj<=0.05,]
DEG_significant = na.omit(DEG_significant)
DEG = na.omit(DEG)
write.csv(DEG, paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Requests/", result_directory, "/DEG.csv"))
write.csv(DEG_significant, paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Requests/", result_directory, "/DEG_significant.csv"))

#Vst normalization and dimensionality reduction analysis
vsd <- vst(dds, blind=FALSE)
pdf(paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/DESeq2/Requests/", result_directory, "/PCA.pdf"), width=5, height=5)
plotPCA(vsd, intgroup=c("R_loops_expression")) + 
  stat_ellipse() + 
  ylim(-200,200) +  
  xlim(-200,200)
plotPCA(vsd, intgroup=c("IDs")) + 
  ylim(-200,200) +  
  xlim(-200,200)
dev.off()