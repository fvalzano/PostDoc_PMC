library(DESeq2)
library(purrr)
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

RNAseq_merge = reduce(RNAseq_runs, function(x, y) merge(x, y, by = "Gene_name"))

