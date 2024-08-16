#Execute in Jobhopper
.libPaths(.libPaths()[3])
library(infercnv)
library(Seurat)

scrna_mimic = readRDS("/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/scrna_mimic_all.rds")
counts_matrix = GetAssayData(scrna_mimic, layer="counts")
barcodes = colnames(scrna_mimic)
annotation = as.data.frame(scrna_mimic$SCT_snn_res.0.4)
annotation = as.matrix(annotation)
infercnv_objs = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                          annotations_file=annotation,
                                          delim="\t",
                                          gene_order_file="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/hg38_gencode_v27.txt",
                                          ref_group_names=c("10", "21", "12", "23", "30", "24"))
infercnv_mimic = infercnv::run(infercnv_objs,
                               cutoff=0.1, #Specified for 10X genomics
                               out_dir="/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCelltype", 
                               cluster_by_groups=TRUE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               resume_mode = T)
scrna_mimic = infercnv::add_to_seurat(seurat_obj = scrna_mimic,
                                           infercnv_output_path = "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCelltype",
                                           top_n = 15)
write_rds(scrna_mimic, "/hpc/pmc_kool/fvalzano/Rstudio_Test1/MIMIC/Data/All/infercnv/CNV_output_xCelltype/scrna_mimic_infercnv.rds")
