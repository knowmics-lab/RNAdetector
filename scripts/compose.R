#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tximport)
  library(GenomicRanges)
  library(rtracklayer)
})

map.ids  <- function (tbl, mapfile) {
  if (!is.na(mapfile) && file.exists(mapfile)) {
    map.table <- read.table(file = mapfile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(map.table) <- c("circRNA_ID", "mapped_id")
    tbl <- tbl %>% left_join(map.table, by="circRNA_ID")
  } else {
    tbl$mapped_id <- NA
  }
  return (tbl)
}

build.granges <- function (file, bedfile, mapfile) {
  m <- read.table(file, stringsAsFactors = FALSE, sep = "\t", skip = 1)
  m <- m[,-12]
  colnames(m) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads","SM_MS_SMS",
                   "non_junction_reads","junction_reads_ratio","circRNA_type","gene_id","strand")
  m <- map.ids(m, mapfile)
  m <- m[order(m$junction_reads,decreasing=T),]
  g <- GRanges(seqnames = Rle(m$chr), ranges = IRanges(m$circRNA_start,m$circRNA_end),strand=m$strand)
  g$id        <- m$circRNA_ID
  g$juncread  <- m$junction_reads
  g$gene      <- m$gene_id
  g$mapped_id <- m$mapped_id
  if (!is.na(bedfile) && file.exists(bedfile)) {
    gr_obj  <- import(bedfile)
    o <- findOverlaps(g, gr_obj, type = "equal")
    g$id[queryHits(o)] <- gr_obj$name[subjectHits(o)]
  }
  return (g)
}

merge.ciri <- function (files, bedfiles, mapfiles) {
  nfile <- length(files)
  glist <- setNames(lapply(1:nfile, function(i) (build.granges(files[i], bedfiles[i], mapfiles[i]))), names(files))
  glist <- GRangesList(glist)
  gcirc <- unlist(glist)
  gcirc$juncread <- NULL
  gcirc.uniq <- unique(gcirc)
  ngcirc <- length(gcirc.uniq)
  gcirc <- gcirc.uniq
  juncread <- matrix(0, nrow = ngcirc, ncol = nfile, dimnames = list(NULL, names(files)))
  for(ifile in 1:nfile){
    tmp <- BiocGenerics::match(glist[[ifile]],gcirc)
    juncread[tmp, ifile] <- glist[[ifile]]$juncread
  }
  juncread <- as.data.frame(juncread)
  return (data.frame(
    id=gcirc$id,
    mapped_id=gcirc$mapped_id,
    name=gcirc$id,
    chr=as(seqnames(gcirc), "vector"),
    start=start(ranges(gcirc)),
    end=end(ranges(gcirc)),
    strand=as.vector(strand(gcirc)),
    length=width(ranges(gcirc)),
    juncread
  ))
}

read.input <- function (input.file, raw = FALSE) {
  return (read.table(input.file, header = FALSE, sep="\t", stringsAsFactors = FALSE))
}

input.as.vector <- function (input, column = 2) {
  if (column > ncol(input)) {
    col <- rep(NA, nrow(input))
  } else {
    col <- input[,column]
  }
  return (setNames(col, input[,1]))  
}

read.samples <- function (input) {
  return (setNames(lapply(input, function (x) {
    table <- read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    no.pos = FALSE
    if (!all("start" %in% colnames(table))) {
      table$start <- NA
      no.pos = TRUE
    }
    if (!all("start" %in% colnames(table))) {
      table$end <- NA
      no.pos = TRUE
    }
    if (!all("length" %in% colnames(table))) {
      if (no.pos) {
        table$length <- NA
      } else {
        table$length <- abs(table$end - table$start) + 1 
      }
    }
    return (table)
  }), names(input)))
}

merge.tables <- function (samples) {
  final <- samples[[1]][,c("id", "mapped_id", "name", "chr", "start", "end", "strand", "length")]
  for (i in 1:length(samples)) {
    tmp <- samples[[i]][,c("id", "counts")]
    colnames(tmp) <- c("id", names(samples)[i])
    final <- final %>% full_join(tmp, by="id")
  }
  final[,8:ncol(final)][is.na(final[,8:ncol(final)])] <- 0
  missing <- which(is.na(final$mapped_id) & is.na(final$name) & 
                     is.na(final$chr) & is.na(final$start) & is.na(final$end) & 
                     is.na(final$strand) & is.na(final$length))
  if (length(missing) > 0) {
    fill.cols  <- c("mapped_id", "name", "chr", "start", "end", "strand", "length")
    missing.id <- final$id[missing]
    for (i in 2:length(samples)) {
      if (length(missing) <= 0) break()
      found      <- which(missing.id %in% samples[[i]]$id)
      found.rows <- which(samples[[i]]$id %in% missing.id[found])
      final[missing[found], fill.cols] <- samples[[i]][found.rows, fill.cols]
      missing    <- missing[-found]
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
  make_option(c("-s", "--toutput"), type="character", default=NULL, help="transcripts expression output file", metavar="character")
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

input <- read.input(opt$input)
if (opt$ciri) {
  merged.input <- merge.ciri(input.as.vector(input, 3), input.as.vector(input, 4), input.as.vector(input, 5))
  write.table(merged.input, file = opt$goutput, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
} else {
  merged.input <- merge.tables(read.samples(input.as.vector(input)))
  write.table(merged.input, file = opt$goutput, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
  if (opt$transcripts) {
    merged.input <- merge.tables(read.samples(input.as.vector(input, column = 4)))
    write.table(merged.input, file = opt$toutput, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}







