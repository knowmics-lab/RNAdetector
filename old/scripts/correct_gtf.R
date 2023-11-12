#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

read.gtf <- function (gtf.file, sep = "\t")  {
  gtf.input <- readr::read_delim(gtf.file, delim = sep, 
                                 col_names = FALSE, comment = "#")
  gffNames <- c("seqname", "source", "feature", "start", "end", 
                "score", "strand", "frame", "attribute")
  names(gtf.input)[1:ncol(gtf.input)] <- gffNames[1:ncol(gtf.input)]
  if (ncol(gtf.input) > 9) 
    stop("The gff file format can not store more than 9 columns!", call.=FALSE)
  return(gtf.input)
}


write.gtf <- function (gtf.data, file.name = "file.gtf") {
  write.table(gtf.data, file.name, sep ="\t", col.names = FALSE, 
              row.names = FALSE, quote = FALSE)
}

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input GTF file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output GTF file", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input) || !file.exists(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required!", call.=FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file is required!", call.=FALSE)
}

gtf <- read.gtf(opt$input)
gtf$attribute <- sapply(strsplit(gtf$attribute, split = ";\\s*", perl = TRUE), function (x) (
  paste0(sapply(strsplit(x, "\\s+"), function(y){
    y <- gsub("^\"|\"$", "", y, perl = TRUE)
    y[2] <- paste0("\"", y[2], "\"")
    return (paste0(y, collapse = " "))
  }), collapse = "; ")
))
write.gtf(gtf, opt$output)
