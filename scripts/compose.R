#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(tximport)

read.gtf <- function (input.gtf, level = NULL, filter.features=NULL, one.by.one = FALSE) {
  a        <- read.table(input.gtf,  stringsAsFactors = FALSE, sep = "\t")
  if (!is.null(level)) {
    if (one.by.one) {
      found <- FALSE
      for (l in level) {
        tmp <- a[a[,3] == l,]
        if (nrow(tmp) > 0) {
          found <- TRUE
          a <- tmp
          break
        }
      }
      if (!found) return (NULL)
    } else {
      a <- a[a[,3] %in% level,]
    }
  }
  if (nrow(a) == 0) return (NULL)
  features <- lapply(strsplit(a[,9],";\\s*", perl = TRUE), function (x) {
    tmp <- strsplit(x, "\\s+\"?", perl = TRUE)
    return (setNames(sapply(tmp, function (x) (ifelse(length(x) > 1, x[2], NA))), sapply(tmp, function (x) (x[1]))))
  })
  all.features <- names(features[[1]])
  if (!is.null(filter.features)) {
    all.features <- intersect(all.features, filter.features)
  }
  columns <- setNames(lapply(all.features, function (f) (sapply(features, function (x) (unname(x[f]))))), all.features)
  a <- a[,-c(2,6,8,9)]
  colnames(a) <- c("chr", "level", "start", "end", "strand")
  for (c in names(columns)) {
    if (!(c %in% colnames(a))) {
      a[[c]] <- columns[[c]]
    }
  } 
  for (c in setdiff(filter.features, all.features)) {
    if (!(c %in% colnames(a))) {
      a[[c]] <- NA
    }
  }
  return (a)
}

harmonize.ciri <- function (input.file, output.file) {
  m <- read.table(input.file, stringsAsFactors = FALSE, sep = "\t", skip = 1)
  m <- m[,c(1,1,2,3,4,11,5)]
  colnames(m) <- c("id", "name", "chr", "start", "end", "strand", "counts")
  m$length <- (m$end - m$start) + 1
  m <- m[,c("id", "name", "chr", "start", "end", "strand", "length", "counts")]
  write.table(m, file = output.file, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
}

read.input <- function (input.file, raw = FALSE) {
  return (read.table(input.file, header = FALSE, sep="\t", stringsAsFactors = FALSE))
}

input.as.vector <- function (input, column = 2) {
  return (setNames(input[,column], input[,1]))  
}

read.samples <- function (input) {
  return (setNames(lapply(input, function (x) (read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t"))), names(input)))
}

merge.tables <- function (samples) {
  final <- samples[[1]][,c("id", "name", "chr", "start", "end", "strand", "length")]
  for (i in 1:length(samples)) {
    tmp <- samples[[i]][,c("id", "counts")]
    colnames(tmp) <- c("id", names(samples)[i])
    final <- final %>% full_join(tmp, by="id")
  }
  missing <- which(is.na(final$name) & is.na(final$chr) & is.na(final$start) & is.na(final$end) & is.na(final$strand) & is.na(final$length))
  if (length(missing) > 0) {
    missing.id <- final$id[missing]
    for (i in 2:length(samples)) {
      if (length(missing) <= 0) break()
      found      <- which(missing.id %in% samples[[i]]$id)
      found.rows <- which(samples[[i]]$id %in% missing.id[found])
      final[missing[found], c("name", "chr", "start", "end", "strand", "length")] <- samples[[i]][found.rows, c("name", "chr", "start", "end", "strand", "length")]
      missing <- missing[-found]
      missing.id <- missing.id[-found]
    }
  }
  return (final)
}

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file", metavar="character"),
  make_option(c("-c", "--ciri"), type="logical", default=FALSE, help="process CIRI output", action = "store_true"),
  make_option(c("-t", "--transcripts"), type="logical", default=FALSE, help="process transcripts output", action = "store_true"),
  make_option(c("-o", "--goutput"), type="character", default=NULL, help="gene expression output file", metavar="character"),
  make_option(c("-t", "--toutput"), type="character", default=NULL, help="transcripts expression output file", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input) || !file.exists(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required!", call.=FALSE)
}

if (is.null(opt$goutput)) {
  print_help(opt_parser)
  stop("Output file is required!", call.=FALSE)
}

if (opt$transcripts && is.null(opt$toutput)) {
  print_help(opt_parser)
  stop("Transcripts output file is required!", call.=FALSE)
}

if (opt$ciri) {
  
} else {
  input <- read.input(opt$input)
  merged.input <- merge.tables(read.samples(input.as.vector(input)))
  write.table(merged.input, file = opt$goutput, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
  if (opt$transcripts) {
    merged.input <- merge.tables(read.samples(input.as.vector(input, column = 3)))
    write.table(merged.input, file = opt$toutput, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}







