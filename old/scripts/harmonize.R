#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tximport)
  library(GenomicRanges)
  library(rtracklayer)
})

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

map.ids  <- function (tbl, map = NULL) {
  if (is.null(map)) {
    mapped           <- tbl
    mapped$mapped_id <- NA
  } else {
    mapped <- tbl %>% left_join(map, by="id")
  }
  mapped <- mapped[,c("id", "mapped_id", "name", "chr", "start", "end", "strand", "length", "counts")]
  return(mapped)
}

harmonize.ciri <- function (input.file, bedfile, output.file, map.table) {
  m <- read.table(input.file, stringsAsFactors = FALSE, sep = "\t", skip = 1)
  m <- m[,c(1,1,2,3,4,11,5)]
  colnames(m) <- c("id", "name", "chr", "start", "end", "strand", "counts")
  m$length <- NA
  m <- m[,c("id", "name", "chr", "start", "end", "strand", "length", "counts")]
  if (!is.null(bedfile) && file.exists(bedfile)) {
    g <- GRanges(seqnames = Rle(m$chr), ranges = IRanges(m$start,m$end),strand=m$strand)
    g$id       <- m$id
    g$gene <- m$gene_id
    gr_obj  <- import(bedfile)
    o <- findOverlaps(g, gr_obj, type = "equal")
    m$id[queryHits(o)] <- gr_obj$name[subjectHits(o)]
  }
  m <- map.ids(m, map.table)
  write.table(m, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
}

harmonize.ciri.quant <- function (input.file, output.file, map.table) {
  m <- read.gtf(input.file)
  m$length <- (m$end - m$start) + 1
  m <- m[,c("circ_id", "circ_id", "chr", "start", "end", "strand", "length", "bsj")]
  colnames(m) <- c("id", "name", "chr", "start", "end", "strand", "length", "counts")
  m <- map.ids(m, map.table)
  write.table(m, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
}

harmonize.htseq <- function (input.file, input.gtf, output.file, map.table) {
  m <- read.table(input.file, stringsAsFactors = FALSE, sep = "\t")
  g <- read.gtf(input.gtf, level=c("gene", "transcript", "exon"), filter.features = c("gene_id", "gene_name"), one.by.one = TRUE)
  colnames(m) <- c("id", "counts")
  m <- m[!grepl(pattern = "^__", m$id, perl = TRUE),]
  if (is.null(g)) {
    m$name <- m$id
    m$chr <- NA
    m$start <- NA
    m$end <- NA
    m$strand <- NA
  } else {
    m <- m %>% left_join(g, by = c("id" = "gene_id"))
    m$name <- m$gene_name
  }
  m$length <- (m$end - m$start) + 1
  m <- m[,c("id", "name", "chr", "start", "end", "strand", "length", "counts")]
  m <- map.ids(m, map.table)
  write.table(m, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
}

harmonize.featureCounts <- function (input.file, input.gtf, output.file, map.table) {
  m <- read.table(input.file, stringsAsFactors = FALSE, sep = "\t")
  m <- m[-1,]
  m[,7] <- as.numeric(m[,7])
  m <- m[,c(1,7)]
  colnames(m) <- c("id", "counts")
  g <- read.gtf(input.gtf, level=c("gene", "transcript", "exon"), filter.features = c("gene_id", "gene_name"), one.by.one = TRUE)
  if (is.null(g)){
    m$name <- m$id
    m$chr <- NA
    m$start <- NA
    m$end <- NA
    m$strand <- NA
  } else {
    m <- m %>% left_join(g, by = c("id" = "gene_id"))
    m$name <- m$gene_name
  }
  m$length <- (m$end - m$start) + 1
  m <- m[,c("id", "name", "chr", "start", "end", "strand", "length", "counts")]
  m <- map.ids(m, map.table)
  write.table(m, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
}

harmonize.salmon <- function (input.file, input.annotation = NULL, output.file, output.transcripts, map.table) {
  m <- read.table(input.file, stringsAsFactors = FALSE, sep = "\t")
  ids <- strsplit(m[-1,1], "|", fixed = TRUE)
  val <- function (x,i) (ifelse(length(x) >= i, x[i], NA))
  annotation <- data.frame(
    seq_id=m[-1,1],
    tx_id=sapply(ids, function(x)(val(x,1))),
    gene_id=sapply(ids, function(x)(val(x,2))),
    tx_name=sapply(ids, function(x)(val(x,5))),
    gene_name=sapply(ids, function(x)(val(x,6)))
  )
  has_genes <- !all(is.na(annotation$gene_id))
  tx2tx <- data.frame(txid=annotation$seq_id,geneid=annotation$tx_id)
  if (has_genes) {
    tx2gene <- data.frame(txid=annotation$seq_id,geneid=annotation$gene_id)
  }
  imported_tx <- tximport(input.file, type = "salmon", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2tx)
  if (has_genes) {
    imported_gn <- tximport(input.file, type = "salmon", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene)
  }
  df_tx <- data.frame(
    id=rownames(imported_tx$counts),
    counts=round(imported_tx$counts),
    length=round(imported_tx$length[rownames(imported_tx$counts),1]),
    chr=NA,
    start=NA,
    end=NA,
    strand=NA,
    row.names = NULL
  ) %>% left_join(annotation, by=c("id"="tx_id"))
  df_tx <- unique(df_tx[,c("id", "tx_name", "chr", "start", "end", "strand", "length", "counts")])
  colnames(df_tx) <- c("id", "name", "chr", "start", "end", "strand", "length", "counts")
  df_tx <- map.ids(df_tx, map.table)
  write.table(df_tx, file = output.transcripts, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)
  if (has_genes) {
    df_gn <- data.frame(
      id=rownames(imported_gn$counts),
      counts=round(imported_gn$counts),
      length=round(imported_gn$length[rownames(imported_gn$counts),1]),
      chr=NA,
      start=NA,
      end=NA,
      strand=NA,
      row.names = NULL
    ) %>% left_join(annotation, by=c("id"="gene_id"))
    df_gn <- unique(df_gn[,c("id", "gene_name", "chr", "start", "end", "strand", "length", "counts")])
    colnames(df_gn) <- c("id", "name", "chr", "start", "end", "strand", "length", "counts")
    df_gn <- map.ids(df_gn, map.table)
    write.table(df_gn, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)   
  }
}

harmonize.stringtie <- function (input.gtf, gc.file, tc.file, output.file, output.transcripts, map.table) {
  g  <- read.gtf(input.gtf, level=c("transcript", "exon"), filter.features = c("gene_id", "transcript_id", "ref_gene_name"), one.by.one = TRUE)
  gc <- read.csv(gc.file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  colnames(gc) <- c("id", "counts")
  tc <- read.csv(tc.file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  colnames(tc) <- c("id", "counts")
  gc$id <- sapply(strsplit(gc$id, split = "|", fixed = TRUE), function (x)(x[1]))
  gc <- gc %>% left_join(g, by=c("id"="gene_id"))
  gc$length <- gc$end - gc$start + 1
  gc <- unique(gc[,c("id", "ref_gene_name", "chr", "start", "end", "strand", "length", "counts")])
  colnames(gc) <- c("id", "name", "chr", "start", "end", "strand", "length", "counts")
  tc <- tc %>% left_join(g, by=c("id"="transcript_id"))
  tc$length <- tc$end - tc$start + 1
  tc <- unique(tc[,c("id", "ref_gene_name", "chr", "start", "end", "strand", "length", "counts")])
  colnames(tc) <- c("id", "name", "chr", "start", "end", "strand", "length", "counts")
  gc <- map.ids(gc, map.table)
  tc <- map.ids(tc, map.table)
  write.table(gc, file = output.file, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)   
  write.table(tc, file = output.transcripts, quote = FALSE, sep = "\t", append = FALSE, row.names = FALSE, col.names = TRUE)   
}

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file", metavar="character"),
  make_option(c("-g", "--annotation"), type="character", default=NULL, help="input gtf file", metavar="character"),
  make_option(c("-a", "--algorithm"), type="character", default=NULL, help="counting algorithm", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL, help="entrez id map", metavar = "character"),
  make_option(c("-o", "--goutput"), type="character", default=NULL, help="gene expression output file", metavar="character"),
  make_option(c("-t", "--toutput"), type="character", default=NULL, help="transcripts expression output file", metavar="character"),
  make_option(c("--stringtie-gene-counts"), type="character", default=NULL, help="stringtie gene counts file created by prepDE.py", metavar="character"),
  make_option(c("--stringtie-tx-counts"), type="character", default=NULL, help="stringtie transcripts counts file created by prepDE.py", metavar="character")
); 

valid.algos   <- c("ciri", "ciriquant", "htseq", "featurecounts", "salmon", "stringtie")
require.annot <- c("htseq", "featurecounts")
require.txout <- c("salmon", "stringtie")

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (opt$help) {
  print_help(opt_parser)
}

if (is.null(opt$input) || !file.exists(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required!", call.=FALSE)
}

if (is.null(opt$algorithm) || !(opt$algorithm %in% valid.algos)) {
  print_help(opt_parser)
  stop("Algorithm is not valid!", call.=FALSE)
}

if ((is.null(opt$annotation) || !file.exists(opt$annotation)) && (opt$algorithm %in% require.annot)) {
  print_help(opt_parser)
  stop("The algorithm requires a valid annotation file!", call.=FALSE)
}

if (is.null(opt$goutput)) {
  print_help(opt_parser)
  stop("The algorithm requires an output file!", call.=FALSE)
}

if (is.null(opt$toutput) && (opt$algorithm %in% require.txout)) {
  print_help(opt_parser)
  stop("The algorithm requires a transcripts output file!", call.=FALSE)
}

if (opt$algorithm == "stringtie" && (is.null(opt[["stringtie-gene-counts"]]) || !file.exists(opt[["stringtie-gene-counts"]]))) {
  print_help(opt_parser)
  stop("Gene counts file is required for stringtie!", call.=FALSE)
}

if (opt$algorithm == "stringtie" && (is.null(opt[["stringtie-tx-counts"]]) || !file.exists(opt[["stringtie-tx-counts"]]))) {
  print_help(opt_parser)
  stop("Transcripts counts file is required for stringtie!", call.=FALSE)
}

map.table <- NULL
if (!is.null(opt$map) && opt$map != "" && file.exists(opt$map)) {
  map.table <- read.table(file = opt$map, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(map.table) <- c("id", "mapped_id")
}

if (opt$algorithm == "ciri") {
  harmonize.ciri(input.file = opt$input, bedfile = opt$annotation, 
                 output.file = opt$goutput, map.table = map.table)
} else if (opt$algorithm == "ciriquant") {
  harmonize.ciri.quant(input.file = opt$input, output.file = opt$goutput, 
                       map.table = map.table)
} else if (opt$algorithm == "htseq") {
  harmonize.htseq(input.file = opt$input, input.gtf = opt$annotation, 
                  output.file = opt$goutput, map.table = map.table)
} else if (opt$algorithm == "featurecounts") {
  harmonize.featureCounts(input.file = opt$input, input.gtf = opt$annotation, 
                          output.file = opt$goutput, map.table = map.table)
} else if (opt$algorithm == "salmon") {
  harmonize.salmon(input.file = opt$input, output.file = opt$goutput, 
                   output.transcripts = opt$toutput, map.table = map.table)
} else if (opt$algorithm == "stringtie") {
  harmonize.stringtie(input.gtf = opt$input, gc.file = opt[["stringtie-gene-counts"]], 
                      tc.file = opt[["stringtie-tx-counts"]], 
                      output.file = opt$goutput, 
                      output.transcripts = opt$toutput, map.table = map.table)
} else {
  print_help(opt_parser)
  stop(paste0("Algorithm ", opt$algorithm, " is not valid!"), call.=FALSE)
}
