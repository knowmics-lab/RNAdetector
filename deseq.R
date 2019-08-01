#### Differential Expression analysis - DESeq 2 ####
#
# - Formatted raw counts input files path. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
# - Output path. Es. "D:/NMSpipeline/DE_results/deseq_test.txt"
# - Path sample info matrix. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/sample_info.txt"
# 
####################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(readr)


deseq_DE_analysis <- function(path_input, path_output, sample_info_path){
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
  
  DESeq.ds<-DESeqDataSetFromMatrix(countData = merged_raw_counts_data,
                                   colData = sample_info,
                                   design = ~ condition)
  
  DESeq.ds <- estimateSizeFactors(DESeq.ds)
  counts_normalized <- counts(DESeq.ds, normalized=TRUE)
  log.norm.counts <- log2(counts_normalized + 1)
  
  DESeq.ds.rlog <- rlog(DESeq.ds, blind = TRUE)
  rlog.norm.counts <- assay(DESeq.ds.rlog)
  
  colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "ctrl")
  DESeq.ds <- DESeq(DESeq.ds)
  DGE.results <- results(DESeq.ds,independentFiltering = TRUE, alpha = 0.05)
  DGE.results.sorted <- DGE.results[order(DGE.results$padj),]
  final_results <- as.data.frame(DGE.results.sorted@listData)
  row.names(final_results) <- DGE.results.sorted@rownames
  final_results <- final_results[,c(2,5,6)]
  write.table(final_results,
              file = path_output,
              row.names = TRUE,
              col.names = TRUE,
              quote = FALSE,
              sep = "\t")
  
}
