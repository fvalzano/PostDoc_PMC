library(DESeq2)
library(pheatmap)
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
Plotting_directory = "HD_Poster/"
Bulk_RNA_files=list()
Bulk_RNA_files = list.files("/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests/20240829_Chris")
#Read the bulkRNA seq runs
Bulk_RNA = list()
for(i in Bulk_RNA_files) {
    Bulk_RNA[[i]] = read.table(paste0("/hpc/pmc_kool/fvalzano/pipelines_fv_output/Data_fetching/Requests/20240829_Chris/", i))
    Bulk_RNA[[i]]= Bulk_RNA[[i]][, c("V2", "V11")]
    colnames(Bulk_RNA[[i]]) = c("Counts", "Gene_name")
}
#Something is weird with the first run - probably someone modified it and saved it?
#Remove the extra row from first BulkRNA run
Bulk_RNA[[1]][1,] = NA
Bulk_RNA[[1]] = na.omit(Bulk_RNA[[1]])
#Check that number of rows now is equal to the rest of bulkRNA runs (60357)
nrow(Bulk_RNA[[1]])
Bulk_RNA_merged = list()
Bulk_RNA_merged = Bulk_RNA[[1]]
for (i in 2:length(Bulk_RNA)) {
  Bulk_RNA_merged <- cbind(Bulk_RNA_merged, Bulk_RNA[[i]])
}
##Remove Gene_name columns duplicates
Bulk_RNA_merged = as.data.frame(Bulk_RNA_merged, row.names = Bulk_RNA_merged$Gene_name)
Bulk_RNA_merged = Bulk_RNA_merged[, -grep("Gene_name", names(Bulk_RNA_merged))]

#Retrieve MB subgroup information - Info was provided by Chris
PMCID=c("PMCID098AAT","PMCID198AAQ","PMCID222AAS","PMCID251AAL","PMCID312AAS","PMCID417AAA",'PMCID432AAR','PMCID481AAQ','PMCID612AAD','PMCID744AAS','PMCID745AAS','PMCID750AAO','PMCID773AAJ','PMCID806AAS','PMCID938AAQ','PMCID961AAR',	'PMCID973AAJ','PMCID981AAH','PMCID053AAT','PMCID072AAQ','PMCID125AAS','PMCID152AAR','PMCID163AAR',	'PMCID170AAS','PMCID179AAB','PMCID190AAQ','PMCID202AAA','PMCID221AAM','PMCID324AAP','PMCID327AAQ','PMCID390AAS','PMCID403AAR','PMCID451AAA','PMCID461AAQ','PMCID495AAO','PMCID592AAS','PMCID664AAR','PMCID739AAO','PMCID778AAB','PMCID781AAS','PMCID816AAI','PMCID879AAO','PMCID896AAJ','PMCID973AAQ','PMCID041AAL','PMCID059AAS','PMCID090AAM','PMCID211AAT','PMCID289AAS','PMCID445AAQ','PMCID524AAQ','PMCID535AAQ','PMCID602AAM','PMCID813AAO','PMCID913AAE','PMCID935AAR','PMCID943AAS','PMCID079AAK','PMCID136AAM','PMCID239AAA','PMCID453AAR','PMCID516AAO','PMCID622AAR','PMCID738AAS','PMCID894AAJ','PMCID926AAR','PMCID457AAT')
##Delete prefix
PMCID = gsub("PMCID", "", PMCID)
Subgroup = c('G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G3','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','G4','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','SHH','Wnt','Wnt','Wnt','Wnt','Wnt','Wnt','Wnt','Wnt','Wnt','Wnt')
PMCID_Subgroup = data.frame(PMCID, Subgroup)

#Load general overview to get the Biomaterial ID matching the PMCID
#One ID has a double entry, delete for now the entry with no info and use the one we are sure is a primary tumor
Overview = read_xlsx("/hpc/pmc_kool/fvalzano/Overview/General_overview_patient_seqdata_ForHPC.xlsx")
Overview=Overview %>% 
  mutate(across(where(is.character), str_trim))
print(Overview[,1], n = 70)
#Overview[c(30,29),] = NA
PMCID_Subgroup_merged = merge(PMCID_Subgroup, Overview, by = "PMCID")
PMCID_Subgroup_merged = PMCID_Subgroup_merged[, c("PMCID", "Subgroup", "RNAseq Patient Biomaterial")]
PMCID_Subgroup_merged = na.omit(PMCID_Subgroup_merged)
PMCID_Subgroup_merged = unique(PMCID_Subgroup_merged)

#Get the biomaterial ID in the same order as they were fetched from the bioinformatic server
Bulk_RNA_files_split = strsplit(Bulk_RNA_files, "_")
BiomaterialID_ordered = list()
for(i in seq_along(Bulk_RNA_files_split)){
    BiomaterialID_ordered[i] = Bulk_RNA_files_split[[i]][1]
    BiomaterialID_ordered = as.character(BiomaterialID_ordered)
}

#Reorder the working dataframe according to the order of fetching from the bioinformatic server
PMCID_Subgroup_ordered <- PMCID_Subgroup_merged[match(BiomaterialID_ordered, PMCID_Subgroup_merged$'RNAseq Patient Biomaterial'), ]
#Apply names on Bulk_RNA_merged dataframe
colnames(Bulk_RNA_merged) = PMCID_Subgroup_ordered$PMCID
#set up deseq2 design
design = as.data.frame(PMCID_Subgroup_ordered$Subgroup)
colnames(design) = "Subgroup"
rownames(design) = colnames(Bulk_RNA_merged)
#Generate dds object
#First column is somehow converted in chr, convert in integer
Bulk_RNA_merged$'451AAA' = as.integer(Bulk_RNA_merged$'451AAA')
dds <- DESeqDataSetFromMatrix(countData = Bulk_RNA_merged,
                              colData = design,
                              design = ~ Subgroup)

##QC the dds object
smallestGroupSize <- 3
retain <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[retain,]
#Normalize counts with vst
dds_vst = vst(dds, blind=FALSE)
#As the triple combination displayed interesting enrichment results, we will focus on the
#upregulated genes in this combination and screen how they look in the patients.
#For the purpose of the poster we will focus on the subset of genes highlighted in the enrichment results

MDG_vs_MD = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_GLI2_vs_MYCN.DNTP53.csv")
MDG_pattern = MDG_vs_MD[order(MDG_vs_MD$log2FoldChange, decreasing = T),]
MDG_pattern_top = MDG_pattern
#From the genes popping out from the enrichment analysis, only PTCH1, EGFR and GLI2 are significant in the comparison between MDG vs E
dds_vst_subset = dds_vst[rownames(dds_vst) %in% c("PTCH1", "EGFR", "GLI2"),]
Bulk_RNA_vst = as.data.frame(assay(dds_vst_subset))
Bulk_RNA_vst = as.data.frame(Bulk_RNA_vst)

annotation = PMCID_Subgroup_ordered
annotation$'RNAseq Patient Biomaterial' = NULL
rownames(annotation) = annotation$PMCID
annotation$PMCID = NULL
pdf(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/Plots/", Plotting_directory, "Overlay_patients_Heatmap_sub.pdf"), height = 2.5, width = 7.5)
pheatmap(Bulk_RNA_vst, scale = "none", annotation = annotation, cluster_row = T, cluster_col = T, show_colnames = F, color=colorRampPalette(c("white", "pink", "red"))(100))
dev.off()


#To broaden our overview, we need to focus not only on the subset of genes
#coming from the enrichment analysis. This was a good first hint but we need
#to broaden the analysis. In order to do so, we concatenated the two DEG analysis
#performed on MDG. Genes present in MDGvsMD will give us GOI different in the two models
#MDGvsE will give us genes differentially expressed from the "healthy" ctrl.
#Load the two DEG analysis
MDG_vs_MD = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_GLI2_vs_MYCN.DNTP53.csv")
MDG_vs_E = read.csv2("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/MYCN_DNTP53_GLI2_vs_EMPTY.csv")
#Combine the two DEGs to have a combo DEG with only genes present in the two
MDG_combo_genes = MDG_vs_E[MDG_vs_E$X %in% MDG_vs_MD$X,]
#Select upregulated genes in MDG
MDG_pattern = MDG_combo_genes[order(MDG_combo_genes$log2FoldChange, decreasing = T),]
MDG_pattern_onlypos = MDG_pattern[MDG_pattern$log2FoldChange>0,]$X
#Plot heatmap of the genes upregulated on patients data
dds_vst_subset = dds_vst[rownames(dds_vst) %in% MDG_pattern_onlypos,]
Bulk_RNA_vst = as.data.frame(assay(dds_vst_subset))
Bulk_RNA_vst = as.data.frame(Bulk_RNA_vst)
pdf(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/Cerebellum_Development/DESEQ2_Analysis/DESeq2_MYCN_MYCN-DNTP53_MYCN-DNTP53-GLI2/Plots/", Plotting_directory, "Overlay_patients_Heatmap_MDG_unique_upregulated.pdf"), height = 20, width = 12.5)
pheatmap(Bulk_RNA_vst, scale = "row", annotation = annotation, cluster_row = T, cluster_col = T, show_colnames = F, color=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()
