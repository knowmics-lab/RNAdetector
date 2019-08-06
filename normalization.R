####################################################### Normalization of raw counts ################################################
#
# - Insert the path of the folder where there are the formatted read counts files. Es: "D:/NMSpipeline/formatted_read_counts/"
# - Insert the path of the output folder. Es: "D:/NMSpipeline/normalized_counts/"
#
####################################################################################################################################
library(readr)

normalization <- function(path_in, path_out){
  input_files_names <- list.files(path_in, pattern = "\\.txt")
  for (input in input_files_names){
    input_path <- paste(path_in, input, sep = "")
    raw_counts <- read_delim(input_path, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
    normalized_counts <- raw_counts
    total_mapped_reads <- sum(normalized_counts$raw_counts)
    normalized_counts$RPM <- 0
    normalized_counts$upper_quartile <- 0
    colnames(normalized_counts) <- c("name", "raw_counts", "RPM", "upper_quartile")
    # RPM
    normalized_counts$RPM <- normalized_counts$raw_counts*(10^6)
    normalized_counts$RPM <- normalized_counts$RPM/total_mapped_reads
    normalized_counts$RPM <- log2(normalized_counts$RPM + 1)
    # Upper quartile
    normalized_counts$upper_quartile <- (normalized_counts$raw_counts * 300)/quantile(normalized_counts$raw_counts[normalized_counts$raw_counts>0],0.75)
    normalized_counts$upper_quartile <- log2(normalized_counts$upper_quartile + 1)
    # Export
    output_path <- paste(path_out, "normalized_", input, sep = "")
    write.table(normalized_counts,
                file = output_path,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }
}

