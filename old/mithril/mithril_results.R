library(readr)
library(pheatmap)
library("RColorBrewer")

mithril_results <- function(path,output_file){
  folder_path <- path
  file_names <- list.files(path, pattern = "metastasis_pathway.txt")
  file_path <- paste(folder_path, file_names[1], sep = "")
  results <- read_delim(file_path, "\t", escape_double = FALSE, trim_ws = TRUE)
  results <- results[results$pValue < 0.05,, drop=FALSE]
  results <- results[,c(2,6)]
  i <- 2
  while(i <= length(file_names)){
    file_path <- paste(folder_path, file_names[i], sep = "")
    df <- read_delim(file_path, "\t", escape_double = FALSE, trim_ws = TRUE)
    df <- df[df$pValue < 0.05,, drop=FALSE]
    df <- df[,c(2,6)]
    results <- merge(x=results, y=df, by.x = "Pathway Name", by.y = "Pathway Name",
                     all.x = TRUE, all.y=TRUE)
    i <- i + 1
  }
  results[is.na(results)] <- 0
  results$`Pathway Name` <- gsub("- Enriched", "", results$`Pathway Name`)
  colnames(results) <- c("pathway_name", file_names)
  write.table(results,
              file = output_file,
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE,
              sep = "\t")
  rownames(results) <- results$pathway_name
  results <- results[,-c(1)]
  breakList <- seq(-20, 20, by= 0.1)
  p <- pheatmap(results,
                color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breakList)),
                cellwidth = 120,
                cellheight = 7.5,
                scale ="none",
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                display_numbers = FALSE)
}



