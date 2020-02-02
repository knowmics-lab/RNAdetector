#!/usr/bin/env Rscript

install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "readr", "tximport", "optparse", "dplyr", "rtracklayer", "plyr"),ask = FALSE)