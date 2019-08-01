#### Differential Expression analysis - LIMMA ####
#
# - Formatted raw counts input files path. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
# - Output path. Es. "D:/NMSpipeline/DE_results/limma_test.txt"
# - Path sample info matrix. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/sample_info.txt"
# 
####################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)
library(readr)

# Parameters
path_input <- "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
path_output <- "D:/NMSpipeline/DE_results/limma_test.txt"
sample_info_path <- "D:/Human/proj_CLL_CD74_case_study/reads_quantification/sample_info.txt"



# Body of the function

# Import and merge raw counts files
raw_counts_file_names <- list.files(path_input, pattern = "\\.txt")
raw_counts_path <- paste(path_input, raw_counts_file_names[1], sep = "")
raw_counts_data <- read_delim(raw_counts_path, "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(raw_counts_data)  <- c("name", raw_counts_file_names[1])
merged_raw_counts_data <- raw_counts_data$name
merged_raw_counts_data <- cbind(merged_raw_counts_data, raw_counts_data[,2])
for (raw_counts_file in raw_counts_file_names[2:length(raw_counts_file_names)]){
  raw_counts_path <- paste(path_input, raw_counts_file, sep = "")
  raw_counts_data <- read_delim(raw_counts_path, "\t", escape_double = FALSE, trim_ws = TRUE)
  colnames(raw_counts_data)  <- c("name", raw_counts_file)
  merged_raw_counts_data <- cbind(merged_raw_counts_data, raw_counts_data[,2])
}
row.names(merged_raw_counts_data) <- merged_raw_counts_data$merged_raw_counts_data
merged_raw_counts_data <- merged_raw_counts_data[,-c(1)]

# Differential expression analysis
sample_info <- read_delim(sample_info_path, "\t", escape_double = FALSE, trim_ws = TRUE)