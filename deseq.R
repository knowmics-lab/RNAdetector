#### Differential Expression analysis - DESeq 2 ####
#
# - Formatted raw counts input files path. Es. "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
# - Output path. Es. "D:/NMSpipeline/DEresults"
# - Condition
# 
#
####################################################

library(BiocInstaller)
biocLite("DESeq2")
library(DESeq2)
library(readr)

# parametri della funzione
path_input <- "D:/Human/proj_CLL_CD74_case_study/reads_quantification/prova/"
path_output <- "D:/NMSpipeline/DEresults"



# funzione
deseq_DE_analysis <- function(path_input, path_output){
  
}

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
sample_info <- data.frame(condition = , row.names = names(merged_raw_counts_data)) # come definire la condition??

DESeq.ds<-DESeqDataSetFromMatrix(countData = merged_raw_counts_data,
                                 colData = sample_info,
                                 design = ~ condition)

DESeq.ds <- estimateSizeFactors(DESeq.ds)
counts_normalized <- counts(DESeq.ds, normalized=TRUE)
log.norm.counts <- log2(counts_normalized + 1)

DESeq.ds.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.ds.rlog)

colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "IDENTIFICATIVO CAMPIONI CONTROLLO NEL DATAFRAME CON LE INFO")
DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds,independentFiltering = TRUE, alpha = 0.05)
DGE.results.sorted <- DGE.results[order(DGE.results$padj),]




# identify genes with the desired adjusted p-value cut-off
DGEgenes<-rownames(subset(DGE.results.sorted, padj<0.05))



