########################################### Calculate Read Per Million Mapped Reads (RPM) ##########################################
#
# - Insert the path of the folder where there are the formatted read counts files. Es: "D:/NMSpipeline/formatted_read_counts/prova/"
# - Insert the path of the output folder. Es: "D:/NMSpipeline/normalized_counts/prova/"
#
####################################################################################################################################

library(readr)

RPM_calculator <- function(path_in, path_out){
  input_files_names <- list.files(path_in, pattern = "\\.txt")
  for (input in input_files_names){
    input_path <- paste(path_in, input, sep = "")
    raw_counts <- read_delim(input_path, "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
    RPM_counts <- raw_counts
    total_mapped_reads <- sum(RPM_counts$raw_counts)
    RPM_counts$RPM <- 0
    colnames(RPM_counts) <- c("name", "raw_counts", "RPM")
    RPM_counts$RPM <- RPM_counts$raw_counts*(10^6)
    RPM_counts$RPM <- RPM_counts$RPM/total_mapped_reads
    RPM_counts$RPM <- log2(RPM_counts$RPM + 1)
    output_path <- paste(path_out, "RPM_", input, sep = "")
    write.table(RPM_counts,
                file = output_path,
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
  }
}


