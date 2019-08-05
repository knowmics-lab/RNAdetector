#### Format MITHrIL input file ####
#
# - DE genes and miRNA input folder. Es. "D:/NMSpipeline/DE_results/"
# - Output path. Es. 
# 
####################################################

library(readr)
library(org.Hs.eg.db)

#Parameters
path_input <- "D:/NMSpipeline/DE_results/"
path_output <- "D:/NMSpipeline/DE_mithril/"

#function
format_input_mithril <- function(path_input, path_output){
  
}

#body of function
input_files <- list.files(path_input, pattern = "\\.txt")
for (input in input_files){
  input_path <- paste(path_input, input, sep = "")
  outputs_path <- paste(path_output, "MITHrIL_", input, sep = "")
  input_matrix <- read_delim(input_path, "\t", escape_double = FALSE, col_names = FALSE,  na = "NA", trim_ws = TRUE, skip = 1)
}


#testing
input_path <- paste(path_input, input_files[1], sep = "")
outputs_path <- paste(path_output, "MITHrIL_", input_files[1], sep = "")
input_matrix <- read_delim(input_path, "\t", escape_double = FALSE, col_names = FALSE,  na = "NA", trim_ws = TRUE, skip = 1)
annotated_matrix <- select(org.Hs.eg.db,
                           keys = input_matrix$X1,
                           keytype = "UCSCKG",
                           columns = "ENTREZID")
