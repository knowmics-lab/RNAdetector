#### Differential Expression analysis - edgeR ####
#
# - Formatted raw counts input files path. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
# - Output path. Es. "D:/NMSpipeline/DE_results/edgeR_test.txt"
# - Path sample info matrix. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/sample_info.txt"
# 
####################################################
library(edgeR)
library(readr)


edgeR_DE_analysis <- function(path_input, path_output, sample_info_path){
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
  sample_types <- factor(sample_info$condition)
  sample_info_edger <- relevel(sample_types, ref = "ctrl")
  design_matrix <- model.matrix (~sample_info_edger)
  colnames(design_matrix) <- levels(sample_types)
  
  edgeR_DGElist <- DGEList(counts = merged_raw_counts_data , group = sample_info_edger)
  edgeR_DGElist<- calcNormFactors(edgeR_DGElist, method = "TMM")
  edgeR_DGElist <- estimateDisp (edgeR_DGElist, design_matrix)
  edger_fit <- glmFit(edgeR_DGElist , design_matrix)
  edger_lrt <- glmLRT(edger_fit)
  
  results <- topTags(edger_lrt, n = Inf, sort.by = "PValue", adjust.method = "BH")
  final_results <- results$table
  final_results <- final_results[,c(1,4,5)]
  final_results <- final_results[order(final_results$logFC, decreasing = TRUE),]
  write.table(final_results,
              file = path_output,
              row.names = TRUE,
              col.names = TRUE,
              quote = FALSE,
              sep = "\t")
}


