# Load the required libraries
library(Seurat)
library(dplyr)
library(biomaRt)

# Define the directory where the count matrices are stored
counts_dir <- "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/count_matrices"
# List all the count matrix files produced by featureCounts
count_files <- list.files(counts_dir, pattern = "*_featureCounts.txt", full.names = TRUE)
#This deletes the summary files from the directory, needs to be executed just once
#summary_files <- list.files(path = counts_dir, pattern = "\\.summary$", full.names = TRUE)
#if (length(summary_files) > 0) {
#  # Delete each ".summary" file
#  file.remove(summary_files)
#  
#  # Print a message indicating successful deletion
#  cat("Deleted", length(summary_files), "'.summary' files from the directory.\n")
#} else {
#  cat("No '.summary' files found in the directory.\n")
#}

# Initialize an empty list to store data
merged_counts <- list()

# Load each count matrix and create a Seurat object
for (file in count_files) {
  # Extract sample name from file name
  sample_name <- basename(file)
  sample_name <- sub("_featureCounts.txt", "", sample_name)
  
  # Read in the count matrix (adjust as necessary depending on file format)
  count_data <- read.table(file, header = TRUE, row.names = 1, check.names = FALSE)
  
  count_data$Gene_ID = rownames(count_data)
  names(count_data)[6] = sample_name
  # Remove any unwanted metadata columns (featureCounts might add extra columns)
  counts <- count_data[,c("Gene_ID", sample_name)]
  # Add the count data to the list
  merged_counts[[sample_name]] <- counts
}
# Merge all the count matrices based on Gene_ID
merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene_ID", all = TRUE), merged_counts)
# Replace NA values with zeros (if any)
merged_df[is.na(merged_df)] <- 0
# Set Gene_ID as rownames and remove the column
rownames(merged_df) <- merged_df$Gene_ID
merged_df$Gene_ID <- NULL
# Convert ENSG gene IDs to gene symbols using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbols for the ENSG IDs in the data
genes_to_convert <- rownames(merged_df)
conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = genes_to_convert,
                    mart = ensembl)

# Merge the conversion table with the counts matrix
# Keep only rows with valid gene symbols
conversion <- conversion %>% filter(hgnc_symbol != "")
merged_df$ensembl_gene_id <- rownames(merged_df)
merged_df <- merge(merged_df, conversion, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# Replace ENSG IDs with gene symbols
merged_df <- merged_df %>%
  dplyr::select(-ensembl_gene_id) %>%
  relocate(hgnc_symbol, .before = everything()) %>%
  rename(Gene_Symbol = hgnc_symbol)
# Remove NAs
merged_df = na.omit(merged_df)
# Set gene symbols as row names
rownames(merged_df) <- make.unique(merged_df$Gene_Symbol, sep="_")
merged_df <- merged_df[, -1]  # Remove the 'Gene_Symbol' column since it's now rownames

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = merged_df, project = "SmartSeq2", min.cells = 3, min.features = 200)
saveRDS(seurat_obj, file = "/hpc/pmc_kool/fvalzano/Ependymoma_Filbin/model_ss2/seurat_merged.rds")

print("Seurat objects have been created and saved.")
