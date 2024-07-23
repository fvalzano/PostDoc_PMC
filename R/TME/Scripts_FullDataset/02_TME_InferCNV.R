#DISCLAIMER_FV: perform inferCNV on Jobhopper, as VSCode instances have troubles with JAGS <--- Fix this someday


#.libPaths(.libPaths()[3]) #Necessary for JobHopper to find the right library directory
library(Seurat)
library(readr)
library(infercnv)

#WithAnnotation
scrna = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_Integration/scrna_harmony.rds")
#Create annotation for each Tumor Cell clusters and TME clusters to be used as experimental and reference group in InferCNV
scrna$Major_classes_FV = as.factor(ifelse(scrna$SCT_snn_res.0.4 %in% "0", "TME_1",
                                          ifelse(scrna$SCT_snn_res.0.4 %in% "28", "TME_2",
                                               ifelse(scrna$SCT_snn_res.0.4 %in% "35", "TME_3",
                                                 ifelse(scrna$SCT_snn_res.0.4 %in% "41", "TME_4",
                                                 paste0("Tumor_Cells_", scrna$SCT_snn_res.0.4)))))

Idents(scrna) = "Dataset"
scrna_list=list()
counts_matrix = list()
barcodes = list()
annotation = list()
infercnv_objs = list()
infercnv = list()

for(i in unique(scrna$Dataset)) {
  scrna_list[[i]] = subset(scrna, idents = i)
  #IMPORTANT:Infercnv can work only on clusters with at least two cells, therefore, we have to delete clusters with less than this cutoff. 
  #This is normally not a problem since you would perform InferCNV without subsetting in single datasets, but here we split according dataset otherwise the process is too computationally demanding
  #First get the ID of the clusters with less than 2 cells
  cluster_counts <- table(scrna_list[[i]]$Major_classes_FV)
  clusters_to_keep <- names(cluster_counts[cluster_counts >= 2])  
  #Then subset the ID from the scrna object
  Idents(scrna_list[[i]]) = "Major_classes_FV"
  scrna_list[[i]] <- subset(scrna_list[[i]], idents = clusters_to_keep)
  #Subsetting retains factors with 0 cells in the levels, manually drop levels with 0 cells
  scrna_list[[i]]$Major_classes_FV = droplevels(scrna_list[[i]]$Major_classes_FV)
  
  counts_matrix[[i]] = GetAssayData(scrna_list[[i]], layer="counts")
  barcodes[[i]] = colnames(scrna_list[[i]])
  annotation[[i]] = as.data.frame(scrna_list[[i]]$Major_classes_FV)
  annotation[[i]] = as.matrix(annotation[[i]])
  infercnv_objs[[i]] = CreateInfercnvObject(raw_counts_matrix=counts_matrix[[i]],
                                            annotations_file=annotation[[i]],
                                            delim="\t",
                                            gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/hg38_gencode_v27.txt",
                                            ref_group_names=unique(scrna_list[[i]]$Major_classes_FV[scrna_list[[i]]$Major_classes_FV %in% c("TME_1", "TME_2", "TME_3", "TME_4")])) # IMPORTANT: We cannot just put TME_1, TME_2... as certain datasets lack some cluster of the TME
  
  #For memory saving purpose
  counts_matrix[[i]] = NULL
  barcodes[[i]] = NULL
  annotation[[i]] = NULL
  scrna_list[[i]] = NULL
  output_dir = dir.create(paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/", i, "_InferCNV"))
  infercnv[[i]] = infercnv::run(infercnv_objs[[i]],
                                cutoff=0.1,
                                out_dir=paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/", i, "_InferCNV"),
                                #analysis_mode = "subclusters",
                                cluster_by_groups=T,
                                denoise=TRUE,
                                HMM=TRUE,
                                resume_mode = T)
  
  #For memory saving purpose
  infercnv[[i]] = NULL
  infercnv_objs[[i]] = NULL
  #scrna[[i]] = infercnv::add_to_seurat(seurat_obj = scrna[[i]],
  #infercnv_output_path = paste0("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/", i,"_InferCNV"),
  #top_n = 15)
}
