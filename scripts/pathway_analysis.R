#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)  
  library(dplyr)
  library(log4r)
  library(rmarkdown)
  library(DT)
  library(ggplot2)
  library(plotly)
  library(d3heatmap)
})

path.env <- new.env(parent = emptyenv())
assign("VERBOSE", NULL, envir = path.env)
assign("LOGGER", NULL, envir = path.env)
assign("MITHRIL.PATH", "/rnadetector/scripts/resources/pathways/", envir = path.env)
assign("PATHWAY.FILE", "/rnadetector/scripts/resources/pathways/pathways.rds", envir = path.env)
assign("TEMPLATE.PATH", "/rnadetector/scripts/resources/pathways/report_template.Rmd", envir = path.env)
assign("FOOTER.PATH", "/rnadetector/scripts/resources/pathways/footer.html", envir = path.env)
assign("CONTRAST.TEMPLATE.PATH", "/rnadetector/scripts/resources/pathways/report_template_contrast.Rmd", envir = path.env)
assign("PATHWAY.TEMPLATE.PATH", "/rnadetector/scripts/resources/pathways/report_template_pathway.Rmd", envir = path.env)

read.contrast.table <- function (contrast, input.path) {
  fn <- paste0(input.path,"/lists/metaseqr_all_out_",contrast,".txt.gz")
  if (!file.exists(fn)) {
    return (NULL)
  }
  return (read.table(fn, sep = "\t", stringsAsFactors = FALSE, header = TRUE))
}

get.p.values <- function (contrast, use.fdr, meta.results) {
  use.metap <- meta.results$complete$params$meta.p != "none" && 
    length(meta.results$complete$params$statistics) > 1
  has.fdr   <- meta.results$complete$params$adjust.method != "none"
  if (use.fdr && has.fdr) {
    if (use.metap) {
      res <- meta.results$complete$meta.fdr[[contrast]]
    } else {
      res <- meta.results$complete$fdr[[contrast]]
    }
  } else {
    if (use.metap) {
      res <- meta.results$complete$meta.p.value[[contrast]]
    } else {
      res <- meta.results$complete$p.value[[contrast]]
    }
  }
  if (is.matrix(res)) {
    res <- setNames(res[,1], rownames(res))
  }
  return (res)
}

get.lfcs     <- function (table) {
  lfc.column <- grep("^log2_normalized_fold_change", colnames(table))
  if (length(lfc.column) != 1) {
    return (NULL)
  }
  res <- table[,c("gene_id", colnames(table)[lfc.column]), drop = FALSE]
  colnames(res) <- c("id", "lfc")
  return(res)
}

get.significant.lfcs <- function (table, contrast, meta.results, source.data, 
                                  p.cut = 0.05, lfc.cut = 0, use.fdr = TRUE) {
  lfcs   <- get.lfcs(table)
  if (is.null(lfcs)) return (NULL)
  pvs    <- get.p.values(contrast, use.fdr, meta.results)
  lfcs$p <- pvs[lfcs$id]
  lfcs   <- lfcs %>% 
    left_join(source.data, by = "id") %>% 
    select("id", "mapped_id", "lfc", "p")
  all.probes <- lfcs$mapped_id[!is.na(lfcs$mapped_id)]
  return (list(
    lfcs=na.omit(
      lfcs %>% 
        filter(p < p.cut) %>%
        filter(abs(lfc) >= lfc.cut)
    ),
    all.probes=all.probes
  ))
}

run.mithril <- function (contrast, lfcs, organism, mithril.path, data.directory) {
  if (!dir.exists(data.directory)) {
    dir.create(data.directory, recursive = TRUE)
  }
  input.file          <- paste0(data.directory, "/mithril_input_", contrast, ".txt")
  output.file.pathway <- paste0(data.directory, "/mithril_output_pathway_", contrast, ".txt")
  output.file.genes   <- paste0(data.directory, "/mithril_output_genes_", contrast, ".txt")
  if (!file.exists(output.file.genes) && !file.exists(output.file.pathway)) {
    lfcs.table          <- lfcs$lfcs[,c("mapped_id", "lfc")]
    zero.probes         <- lfcs$all.probes[!(lfcs$all.probes %in% lfcs.table$mapped_id)]
    zero.probes         <- cbind(zero.probes, 0)
    colnames(zero.probes) <- colnames(lfcs.table)
    lfcs.table          <- rbind(lfcs.table, zero.probes)
    write.table(lfcs.table, file = input.file, append = FALSE, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = FALSE)
    command.line <- paste0("java -jar ", mithril.path, "/MITHrIL2.jar mithril ",
                           "-m -i '", input.file, "' -organism ", organism, 
                           " -o '",output.file.pathway,"' -p '", output.file.genes, 
                           "' -no-complete")
    code <- system(command.line, ignore.stdout = TRUE, ignore.stderr = TRUE, 
                   wait = TRUE)
    if (code != 0) {
      return (NULL)
    }
  }
  pathway.table <- read.delim(output.file.pathway, stringsAsFactors = FALSE)
  colnames(pathway.table)[1] <- "Pathway.Id"
  pathway.table$Pathway.Name <- gsub("\\s+-\\s+Enriched", "", 
                                     pathway.table$Pathway.Name, perl = TRUE)
  genes.table   <- read.delim(output.file.genes, stringsAsFactors = FALSE)
  colnames(genes.table)[1] <- "Pathway.Id"
  genes.table$Pathway.Name <- gsub("\\s+-\\s+Enriched", "", 
                                   genes.table$Pathway.Name, perl = TRUE)
  genes.table.fdr <- tapply(1:nrow(genes.table), genes.table$Pathway.Id, 
                            function (i) (setNames(p.adjust(
                              genes.table$pValue[i], method = "fdr"), 
                              genes.table$Gene.Id[i])))
  return (list(
    pathway.table=pathway.table,
    genes.table=genes.table,
    genes.table.fdr=genes.table.fdr
  ))
}

get.pathway.data <- function (pathway, organism) {
  if (!exists("PATHWAY.DATA", envir = path.env)) {
    assign("PATHWAY.DATA", readRDS(get("PATHWAY.FILE", envir = path.env)), envir = path.env)
  }
  data <- get("PATHWAY.DATA", envir = path.env)
  pathway <- gsub("path:", "", pathway, fixed = TRUE)
  if (!(organism %in% names(data))) {
    return(NULL)
  }
  if (!(pathway %in% names(data[[organism]]))) {
    return(NULL)
  }
  return (data[[organism]][[pathway]])
}

contrasts.to.text <- function (contrast) {
  return (sapply(strsplit(contrast, "_vs_"), function (s) (paste(s, collapse = " vs "))))
}

pathway.heatmap <- function (path.results, p.cut, use.fdr) {
  path.processed <- setNames(lapply(path.results, function (res) {
    pt <- res$pathway.table
    if (use.fdr) {
      pt$p <- pt$Adjusted.pValue
    } else {
      pt$p <- pt$pValue
    }
    return (pt)
  }), names(path.results))
  all.pathways <- unique(unlist(lapply(path.processed, 
                                       function (p) (p$Pathway.Id[p$p < p.cut]))))
  path.filtered <- setNames(lapply(path.processed, function(p) {
    sel <- p[p$Pathway.Id %in% all.pathways,]
    return (setNames(sel$Corrected.Accumulator, sel$Pathway.Name))
  }), names(path.processed))
  path.matrix <- do.call(cbind, path.filtered)
  colnames(path.matrix) <- contrasts.to.text(colnames(path.matrix))
  return (
    d3heatmap(path.matrix, scale="row", colors="RdBu", dendrogram = "row", 
              xaxis_font_size = "8pt", yaxis_font_size = "8pt")
  )
}

get.genes.table <- function (path.res, contrast, pathway, p.cut) {
  table <- path.res[[contrast]]$genes.table %>% 
    filter(Pathway.Id == pathway) %>%
    left_join(source.data, by = c("Gene.Id"="mapped_id")) %>%
    select("Gene.Id", "Gene.Name", "name", "Perturbation", "pValue")
  fdr   <- path.res[[contrast]]$genes.table.fdr[[pathway]]
  table$Gene.Name[!is.na(table$name)] <- table$name[!is.na(table$name)]
  table$name <- NULL
  if (nrow(table) == 0) return (NULL)
  table$FDR <- fdr[table$Gene.Id]
  table$color.acc <- 1
  table$color.acc[table$pValue < p.cut & table$Perturbation < 0] <- 2
  table$color.acc[table$pValue < p.cut & table$Perturbation > 0] <- 3
  table$name.acc <- "Non significant"
  table$name.acc[table$pValue < p.cut & table$Perturbation < 0] <- "Pert < 0"
  table$name.acc[table$pValue < p.cut & table$Perturbation > 0] <- "Pert > 0"
  return (table)
}

disp <- function(...) {
  verbose <- get("VERBOSE",envir=path.env)
  if (!is.null(verbose) && verbose) {
    message("\n",...,appendLF=FALSE)
  }
  logger <- get("LOGGER",envir=path.env)
  levalias <- c("one","two","three","four","five")
  if (!is.null(logger)) {
    switch(levalias[level(logger)],
           one = { debug(logger,paste0(...)) },
           two = { info(logger,gsub("\\n","",paste0(...))) },
           three = { warn(logger,gsub("\\n","",paste0(...))) },
           four = { error(logger,gsub("\\n","",paste0(...))) },
           five = { fatal(logger,gsub("\\n","",paste0(...))) }
    )
  }
}

stopwrap <- function(...,t="fatal") {
  logger <- get("LOGGER",envir=path.env)
  if (!is.null(logger)) {
    if (t=="fatal")
      fatal(logger,gsub("\\n","",paste0(...)))
    else
      error(logger,gsub("\\n","",paste0(...)))
  }
  stop(paste0(...))
}

warnwrap <- function(...,now=FALSE) {
  logger <- get("LOGGER",envir=path.env)
  if (!is.null(logger))
    warn(logger,gsub("\\n","",paste0(...)))
  if (now)
    warning(paste0(...),call.=FALSE,immediate.=TRUE)
  else
    warning(paste0(...),call.=FALSE)
}

elap2human <- function(start.time) {
  start.time <- as.POSIXct(start.time)
  dt <- difftime(Sys.time(),start.time,units="secs")
  ndt <- as.numeric(dt)
  if (ndt<60)
    format(.POSIXct(dt,tz="GMT"),"%S seconds")
  else if (ndt>=60 && ndt<3600)
    format(.POSIXct(dt,tz="GMT"),"%M minutes %S seconds")
  else if (ndt>=3600 && ndt<86400)
    format(.POSIXct(dt,tz="GMT"),"%H hours %M minutes %S seconds")
  else if (ndt>=86400)
    format(.POSIXct(dt,tz="GMT"),"%d days %H hours %M minutes %S seconds")
}

build.report <- function (
  input.directory,
  output.directory,
  degs.p.cut = 0.05,
  degs.lfc.cut = 0,
  degs.use.fdr = TRUE,
  pathway.organism = "hsa",
  pathway.p.cut = 0.05,
  pathway.use.fdr = TRUE,
  verbose = TRUE,
  run.log = TRUE
) {
  output.paths <- list(
    "logs"=paste0(output.directory,"/logs"),
    "data"=paste0(output.directory,"/data"),
    "libs"=paste0(output.directory,"/libs")
  )
  lapply(output.paths, function(d) (dir.create(d, showWarnings = FALSE, recursive = TRUE)))
  assign("VERBOSE", verbose, envir=path.env)
  if (run.log) {
    logger <- create.logger(
      logfile=file.path(output.paths$logs, "/run.log"), level=2, logformat="%d %c %m")
  } else {
    logger <- NULL
  }
  assign("LOGGER", logger, envir=path.env)
  
  TB <- Sys.time()
  disp(strftime(Sys.time()),": Data processing started...\n")
  disp("DEGs analysis directory: ",input.directory)
  data.file <- paste0(input.directory, "/data/analysis.RData")
  if (!file.exists(data.file)) {
    stopwrap("Unable to find analysis input data.")
  }
  load(data.file, envir = globalenv())
  source.data$mapped_id <<- as.character(source.data$mapped_id)
  disp("Contrasts: ",paste(contrasts.list,collapse=", "))
  disp("DEGs Analysis Parameters:")
  if (degs.use.fdr) {
    disp("  FDR p-value threshold: ", degs.p.cut)
  } else {
    disp("  p-value threshold: ", degs.p.cut)
  }
  disp("  Log-Fold-Change threshold: ", degs.lfc.cut)
  disp("Pathway Analysis Parameters:")
  disp("  Organism: ", pathway.organism)
  if (pathway.use.fdr) {
    disp("  FDR p-value threshold: ", pathway.p.cut)
  } else {
    disp("  p-value threshold: ", pathway.p.cut)
  }
  ncon       <- length(contrasts.list)
  con.tables <- setNames(vector("list", ncon), contrasts.list)
  all.lfcs   <- setNames(vector("list", ncon), contrasts.list)
  path.res   <- setNames(vector("list", ncon), contrasts.list)
  for (con in contrasts.list) {
    disp("Analyzing contrast ", con, ":")
    disp("  Reading contrast table")
    con.table <- read.contrast.table(con, input.directory)
    if (is.null(con.table)) {
      disp("  Unable to contrast table.")
      next()
    }
    disp("  Computing significant DEGs: ")
    sign.lfcs <- get.significant.lfcs(con.table, con, meta, source.data, 
                                      degs.p.cut, degs.lfc.cut, degs.use.fdr)
    if (is.null(sign.lfcs) || !is.list(sign.lfcs)) {
      disp("    Unable to compute significant DEGs.")
      next()
    } else if (nrow(sign.lfcs$lfcs) == 0) {
      disp("    No DEGs found.")
      next()
    }
    disp("    Found ", nrow(sign.lfcs$lfcs), " DEGs.")
    disp("  Running MITHrIL analysis:")
    mith.res <- run.mithril(con, sign.lfcs, pathway.organism, 
                            get("MITHRIL.PATH", envir = path.env), 
                            output.paths$data)
    if (is.null(mith.res) || !is.list(mith.res)) {
      disp("    Unable to run MITHrIL analysis.")
      next()
    } else {
      disp("    MITHrIL analysis completed.")
    }
    con.tables[[con]] <- con.table
    all.lfcs[[con]]   <- sign.lfcs
    path.res[[con]]   <- mith.res
  }
  path.res <- Filter(function (x)(!is.null(x)), path.res)
  c.list   <- names(path.res)
  path.with.p <- setNames(lapply(path.res, function (res) {
    pt <- res$pathway.table
    if (pathway.use.fdr) {
      pt$p <- pt$Adjusted.pValue
    } else {
      pt$p <- pt$pValue
    }
    pt$color.acc <- 1
    pt$color.acc[pt$p < pathway.p.cut & pt$Corrected.Accumulator < 0] <- 2
    pt$color.acc[pt$p < pathway.p.cut & pt$Corrected.Accumulator > 0] <- 3
    pt$name.acc <- "Non significant"
    pt$name.acc[pt$p < pathway.p.cut & pt$Corrected.Accumulator < 0] <- "Acc < 0"
    pt$name.acc[pt$p < pathway.p.cut & pt$Corrected.Accumulator > 0] <- "Acc > 0"
    return (pt)
  }), names(path.res))
  path.processed <- setNames(lapply(path.with.p, function (pt) {
    pt$color.acc <- NULL
    pt$name.acc <- NULL
    return (pt %>% filter(p < pathway.p.cut))
  }), names(path.with.p))
  path.tables <- setNames(lapply(path.with.p, function (pt) {
    pt$color.acc <- NULL
    pt$name.acc <- NULL
    pt$Pathway.Id <- gsub("path:", "", pt$Pathway.Id)
    pt$Raw.Accumulator <- NULL
    pt$Impact.Factor   <- NULL
    colnames(pt) <- c("Id", "Pathway Name", "pOra", "pAcc", "Acc", "pValue", "FDR", "p")
    pt <- pt %>% filter(p < pathway.p.cut)
    pt$p <- NULL
    return (pt)
  }), names(path.with.p))
  disp("Creating HTML report.")
  rmarkdown::render(
    input = get("TEMPLATE.PATH", envir = path.env),
    output_file = "index.html",
    output_dir = output.directory,
    output_options = list(self_contained=FALSE, lib_dir=output.paths$libs, includes=list(after_body=get("FOOTER.PATH", envir = path.env))),
    knit_root_dir = output.directory,
    envir = environment(),
    clean = TRUE,
    quiet = TRUE
  )
  analysis.data <- list(
    contrast.tables=con.tables,
    log.fold.changes=all.lfcs,
    pathway.results=path.res,
    contrasts.list=contrasts.list
  )
  save(analysis.data, file = paste0(output.paths$data, "/analysis.RData"))
  disp("\n",strftime(Sys.time()),": Data processing finished!\n")
  exec.time <- elap2human(TB)
  disp("\n","Total processing time: ",exec.time,"\n\n")
  return(environment())
}

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input directory", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output directory", metavar="character"),
  make_option(c("--degs-p"), type="numeric", default=0.05, help="DEGs: p-value threshold", metavar="number"),
  make_option(c("--degs-lfc"), type="numeric", default=0, help="DEGs: log-fold-change threshold", metavar="number"),
  make_option(c("--degs-no-fdr"), type="logical", default=TRUE, help="DEGs: disable FDR threshold", action = "store_false"),
  make_option(c("--path-org"), type="character", default="hsa", help="Pathway Analysis: organism", metavar="number"),
  make_option(c("--path-p"), type="numeric", default=0.05, help="Pathway Analysis: p-value threshold", metavar="number"),
  make_option(c("--path-no-fdr"), type="logical", default=TRUE, help="Pathway Analysis: disable FDR threshold", action = "store_false")
); 

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$input) || !dir.exists(opt$input)) {
  print_help(opt_parser)
  stop("Input directory is required!", call.=FALSE)
}

if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output directory is required!", call.=FALSE)
}

res <- build.report(
  input.directory = opt$input, 
  output.directory = opt$output, 
  degs.p.cut = opt[["degs-p"]],
  degs.lfc.cut = opt[["degs-lfc"]],
  degs.use.fdr = opt[["degs-no-fdr"]],
  pathway.organism = opt[["path-org"]],
  pathway.p.cut = opt[["path-p"]],
  pathway.use.fdr = opt[["path-no-fdr"]]
)
