#.libPaths(.libPaths()[3]) #Necessary for JobHopper to find the right library directory
library(Seurat)
library(readr)
library(infercnv)

#DISCLAIMER_FV: perform inferCNV on Jobhopper, as VSCode instances have troubles with JAGS <--- Fix this someday

#WithAnnotation
scrna = read_rds("/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/scrna_harmony.rds")
Idents(scrna) = "SCT_snn_res.0.4"
counts_matrix = GetAssayData(scrna, layer="counts")
barcodes = colnames(scrna)
annotation = scrna$SCT_snn_res.0.4
annotation = as.matrix(annotation)
infercnv = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/hg38_gencode_v27.txt",
                                    ref_group_names=c("0", "28", "41", "35"))

infercnv = infercnv::run(infercnv,
                              cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/WithAnnotation", 
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T)

scrna = infercnv::add_to_seurat(seurat_obj = scrna,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/WithAnnotation",
                                           top_n = 15)
write_rds(scrna, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_InferCNV/scrna_infercnv_with_annotations.rds")

#Without Annotation
scrna$AllCells = "Cells"
Idents(scrna) = "AllCells"
counts_matrix = GetAssayData(scrna, layer="counts")
barcodes = colnames(scrna)
annotation = scrna$AllCells
annotation = as.matrix(annotation)
infercnv = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotation,
                                    delim="\t",
                                    gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/hg38_gencode_v27.txt",
                                    ref_group_names=c())

infercnv = infercnv::run(infercnv,
                              cutoff=0.1,
                             out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/NoAnnotation", 
                             cluster_by_groups=F, 
                             denoise=TRUE,
                             HMM=TRUE,
                             resume_mode = T,
                             k_obs_groups = 2)

scrna = infercnv::add_to_seurat(seurat_obj = scrna,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/InferCNV/Output/NoAnnotation",
                                           top_n = 15)
write_rds(scrna, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/TME/TME_files_March24/Seurat_subsets/Post_InferCNV/scrna_infercnv_no_annotations.rds")
