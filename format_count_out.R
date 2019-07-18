##### Change CIRI, HTseq-count, FeatureCounts and Salmon output files in the same output file format #####
#
# -Insert the path of the folder where there are the read counts files. Es. "D:/NMSpipeline/read_counts/"
# -Insert the path of the output folder. Es. "D:/NMSpipeline/formatted_read_counts/"
#
##########################################################################################################

library(readr)
format_output_count_files <- function(path, output_path){
  counts_files_names <- list.files(path, pattern = "\\.txt")
  for (counts in counts_files_names){
    counts_path <- paste(path, counts, sep = "")
    output_counts_path <- paste(output_path, "formatted_", counts, sep = "")
    splitted_name <- unlist(strsplit(counts, "_"))
    suff <- splitted_name[2]
    if (suff == "ht.txt"){
      htseq_data <- read_delim(counts_path, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
      colnames(htseq_data) <- c("name", "raw_counts")
      write.table(htseq_data,
                  file = output_counts_path, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    } else if (suff == "fc.txt") {
      featureCounts_data <- read_delim(counts_path, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 1)
      featureCounts_data <- featureCounts_data[,c(1,7)]
      colnames(featureCounts_data) <- c("name", "raw_counts")
      write.table(featureCounts_data,
                  file = output_counts_path, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}

