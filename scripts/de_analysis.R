#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(rjson)
  library(dplyr)
  library(metaseqR)
})

unslash <- function(dirs) (sub("/+$", "", dirs))

option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL, help="config file", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (opt$help) {
  print_help(opt_parser)
}

if (is.null(opt$config) || !file.exists(opt$config)) {
  print_help(opt_parser)
  stop("Config file is required!", call.=FALSE)
}

config <- fromJSON(file = opt$config)

if (is.null(config$description.file) || !file.exists(config$description.file)) {
  stop("Missing data description file.", call. = FALSE)
}
if (is.null(config$data.file) || !file.exists(config$data.file)) {
  stop("Missing data file.", call. = FALSE)
}
if (is.null(config$data.type)) {
  config$data.type <- "gene"
}
if (is.null(config$conditions.variables) || length(config$conditions.variables) <= 0) {
  stop("Missing condition variables.", call. = FALSE)
}
if (is.null(config$contrasts) || length(config$contrasts) <= 0) {
  stop("Missing contrasts.", call. = FALSE)
}
if (is.null(config$output.directory)) {
  stop("Missing output directory", call. = FALSE)
}
config$output.directory <- unslash(config$output.directory)

descriptions <- read.delim(config$description.file, stringsAsFactors=FALSE)
if (!("SampleName" %in% colnames(descriptions))) {
  stop("Invalid description file.", call. = FALSE)
}
data         <- read.delim(config$data.file, stringsAsFactors=FALSE)
variables    <- config$conditions.variables[config$conditions.variables %in% colnames(descriptions)]
if (length(variables) <= 0) {
  stop("No valid variables found.", call. = FALSE)
}
descriptions$condition <- do.call(paste, 
                                  c(
                                    lapply(variables, function (v)(descriptions[[v]])), 
                                    list(sep="_")))
samples.list <- tapply(descriptions$SampleName, descriptions$condition, function(x) (x))
contrasts.list <- sapply(config$contrasts, function(x)(paste(x, collapse = "_vs_")))

data$gc     <- 0
common.cols <- c("id","name","chr","start","end","strand","length","gc")
if (all("mapped_id" %in% colnames(data))) {
  data$mapped_id <- NULL #### TEMPORANEO DA RIMUOVERE PER PATHWAY ANALYSIS
}
other.cols  <- setdiff(colnames(data), common.cols)
data        <- data[,c(common.cols, other.cols)]

if (is.null(config$parameters)) {
  params <- list()
} else {
  params <- config$parameters
}

check.list <- function (what, defaults=NULL) {
  if (is.null(what) || !is.list(what)) {
    res <- list()
  } else {
    res <- what
  }
  if (!is.null(defaults) && is.list(defaults)) {
    res <- modifyList(defaults, res)
  }
  return (res)
}

check.vector <- function (what, default) {
  if (is.null(what) || !is.vector(what)) {
    res <- default
  } else {
    res <- what
  }
  return (res)
}

restrict.cores <- 1 / check.vector(params$num.cores, 1)
restrict.cores <- max(min(restrict.cores, 1), 1 / detectCores())
pcut <- check.vector(params$pcut, 0.05)
pcut <- max(min(pcut, 1), 0)
when.apply.filter <- check.vector(params$when.apply.filter, "prenorm")
norm.algo <- check.vector(params$norm, "edger")
norm.algo.params <- check.list(params$norm.args, get.defaults("normalization", norm.algo))
if (norm.algo == "deseq") {
  if (!is.null(norm.algo.params$locfunc) && is.character(norm.algo.params$locfunc)) {
    if (norm.algo.params$locfunc == 'shorth') {
      norm.algo.params$locfunc <- genefilter::shorth
    } else {
      norm.algo.params$locfunc <- stats::median
    }
  }
}
stats.algo <- check.vector(params$stats, "limma")
stats.algo.params <- check.list(params$stats.args)
stats.algo.params <- setNames(lapply(stats.algo, function (a) {
  tmp <- check.list(stats.algo.params[[a]], get.defaults("statistics", a))
}), stats.algo)
default.filters <- get.defaults("gene.filter", "")
gene.filters <- check.list(params$filters)
for (f in names(gene.filters)) {
  if (is.null(gene.filters[[f]])) {
    gene.filters[[f]] <- NULL 
  } else {
    gene.filters[[f]] <- modifyList(default.filters[[f]], gene.filters[[f]])
  }
}

suppressWarnings({
  meta <- metaseqr(
    counts=data,
    sample.list=samples.list,
    contrast=contrasts.list,
    id.col = 1,
    name.col = 2,
    gc.col = 8,
    annotation = "embedded",
    org = "custom",
    trans.level = config$data.type,
    count.type = "gene",
    when.apply.filter = when.apply.filter,
    normalization = norm.algo,
    norm.args = norm.algo.params,
    statistics = stats.algo,
    stat.args = stats.algo.params,
    adjust.method = check.vector(params$adjust.method, "qvalue"),
    meta.p = check.vector(params$meta.p.method, "simes"),
    gene.filters = gene.filters,
    qc.plots=c(
      "mds", "biodetection", "countsbio", "saturation", "readnoise",
      "filtered", "correl", "pairwise", "boxplot", "lengthbias", 
      "meandiff", "meanvar", "deheatmap", "volcano", #"rnacomp", 
      "biodist","venn"
    ),
    fig.format=check.vector(params$fig.formats, c("png","pdf")),
    export.where=config$output.directory,
    export.what = c(
      "annotation", "p.value", "adj.p.value", "meta.p.value", 
      "adj.meta.p.value", "fold.change", "stats", "counts", 
      "flags"
    ),
    export.scale = c("natural", "log2", "vst"),
    export.values = c("raw", "normalized"),
    export.stats = c("mean", "median", "sd", "mad", "cv", "rcv"),
    export.counts.table = TRUE,
    out.list = TRUE
  )
  save(data, meta, samples.list, contrasts.list, 
       file = paste0(config$output.directory, "/data/analysis.RData"))
})

