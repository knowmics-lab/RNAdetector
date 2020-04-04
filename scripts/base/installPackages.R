#!/usr/bin/env Rscript

install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "readr", "tximport", 
                        "optparse", "dplyr", "rtracklayer", "plyr", 
                        "survcomp", "VennDiagram", "knitr", "zoo",
                        "devtools", "plotly", "d3heatmap", "rmarkdown", 
                        "DT"),ask = FALSE)

remotes::install_github("alaimos/metaseqR", dependencies = TRUE, upgrade = "always")