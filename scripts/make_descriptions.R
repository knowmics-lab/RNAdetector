#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(plyr)
})

read.descriptions <- function (input, samples.list, default.group) {
  return (c(lapply(input, function (x) { 
    tmp <- read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
    colnames(tmp)[1] <- "SampleName"
    if (!("SampleGroup" %in% colnames(tmp))) {
      tmp$SampleGroup <- default.group
    }
    tmp <- tmp[,c("SampleName", "SampleGroup", colnames(tmp)[-which(colnames(tmp) %in% c("SampleName", "SampleGroup"))])]
    return (tmp)
  }), 
  list(
    data.frame(SampleName=samples.list, SampleGroup=default.group, stringsAsFactors = FALSE)
  )))
}

remove.invalid.samples <- function (input, samples.list) {
  return (lapply(input, function (x) {
    r <- which(!(x$SampleName %in% samples.list))
    if (length(r) > 0) {
      x <- x[-r,]
    }
    return (x)
  }))
}

remove.duplicated.samples <- function (input) {
  contained <- c()
  for (i in 1:length(input)) {
    x <- input[[i]]
    dups <- which(x$SampleName %in% contained)
    if (length(dups) > 0) {
      x <- x[-dups,]
    }
    input[[i]] <- x
    contained <- c(contained, x$SampleName)
  }
  return (input)
}

option_list <- list(
  make_option(c("-d", "--descriptions"), type="character", default=NULL, help="text file containing a list of description matrices", metavar="character"),
  make_option(c("-g", "--group"), type="character", default=FALSE, help="text file containing sample group name", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=FALSE, help="a list of valid sample names", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$descriptions) || !file.exists(opt$descriptions)) {
  print_help(opt_parser)
  stop("Descriptions file is required!", call.=FALSE)
}

if (is.null(opt$group)) {
  print_help(opt_parser)
  stop("Default group name is required!", call.=FALSE)
}

if (is.null(opt$samples) || !file.exists(opt$samples)) {
  print_help(opt_parser)
  stop("Samples list file is required!", call.=FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file is required!", call.=FALSE)
}

input         <- readLines(opt$descriptions)
samples       <- readLines(opt$samples)
default.group <- opt$group
descriptions  <- read.descriptions(input, samples, default.group)
final.df      <- rbind.fill(Filter(function (x) (nrow(x) > 0), remove.duplicated.samples(remove.invalid.samples(descriptions, samples))))
write.table(final.df, file = opt$output, append = FALSE, quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)




