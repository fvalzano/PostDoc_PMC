library(DESeq2)
library(purrr)
library(apeglm)
library(ggplot2)
library(umap)
library(uwot)
#Load bulk RNAseq runs from specific requests folder
##Modify request folder
result_directory = "20240701_Francesco"
IDs = c("806AAS", "222AAS", "745AAS")
IDs = IDs[order(IDs, decreasing = F)]
##Set up list to contain the single bulk RNA seq runs, delete unnecessary columns and rename remaining ones (counts and gene name)
RNAseq_files = list.files(paste0("/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/Requests/", result_directory))
RNAseq_runs = list()
for (file in RNAseq_files) {
   RNAseq_runs[[file]]=read.table(paste0("/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/Data_fetching/Requests/", result_directory, "/", file), )
   RNAseq_runs[[file]]= RNAseq_runs[[file]][, c("V2", "V11")]
   colnames(RNAseq_runs[[file]]) = c("Counts", "Gene_name")
}

##Merge the single runs in one
RNAseq_merged = RNAseq_runs[[1]]
for (i in 2:length(RNAseq_runs)) {
  RNAseq_merged <- cbind(RNAseq_merged, RNAseq_runs[[i]])
}

##Remove Gene_name columns duplicates
RNAseq_merged = as.data.frame(RNAseq_merged, row.names = RNAseq_merged$Gene_name)
RNAseq_merged = RNAseq_merged[, -grep("Gene_name", names(RNAseq_merged))]
##Rename dataframe with IDs and material of origin
IDs_patients = paste0(IDs, "_Patients")
IDs_tumoroid = paste0(IDs, "_Tumoroid")
colnames(RNAseq_merged) = c(IDs_patients, IDs_tumoroid)
##Write analysis design and set up dds object
design = as.data.frame(rep(c("Patient", "Tumoroid"), each = 3))
colnames(design) = "material"
design$IDs = rep(IDs,2)
rownames(design) = colnames(RNAseq_merged)
dds <- DESeqDataSetFromMatrix(countData = RNAseq_merged,
                              colData = design,
                              design = ~ material)

##QC the dds object
smallestGroupSize <- 3
retain <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[retain,]

##Differentially Expressed Gene Analysis with apeglm shrinkage of fold changes
resultsNames(dds)
dds <- DESeq(dds)
res <- results(dds)
resLFC <- lfcShrink(dds, coef="material_Tumoroid_vs_Patient", type="apeglm")
DEG = as.data.frame(resLFC)
DEG = DEG[DEG$padj<=0.05,]
DEG = na.omit(DEG)
write.csv(DEG, paste0("/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2/Requests/", result_directory, "/DEG.csv"))

#Vst normalization and dimensionality reduction analysis
vsd <- vst(dds, blind=FALSE)
pdf(paste0("/hpc/pmc_kool/fvalzano/PostDoc_PMC/pipelines/DESeq2/Requests/", result_directory, "/PCA.pdf"), width=5, height=5)
plotPCA(vsd, intgroup=c("material")) + stat_ellipse()
plotPCA(vsd, intgroup=c("material", "IDs"))
dev.off()