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

limma_DE_analysis <- function(path_input, path_output, sample_info_path){
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
  design_matrix <- model.matrix(~0+sample_info$condition)
  sample_types <- levels(factor(sample_info$condition))
  colnames(design_matrix) <- sample_types

  voomTransformed<-voom(merged_raw_counts_data, design_matrix, plot=FALSE)
  voomed.fitted <- lmFit(voomTransformed, design=design_matrix)

  contrasts <- makeContrasts(sample-ctrl, levels = design_matrix)
  fit_2 <- contrasts.fit(voomed.fitted, contrasts)
  fit_2 <-eBayes(fit_2)

  # Format output
  results <- topTable(fit_2, coef=1, number = Inf, adjust.method ="BH", sort.by="logFC")
  results <- results[,c(1,4,5)]
  results <- results[order(results$logFC, decreasing = TRUE),]
  write.table(results,
              file = path_output,
              col.names = TRUE,
              row.names = TRUE,
              quote = FALSE,
              sep = "\t")
}
