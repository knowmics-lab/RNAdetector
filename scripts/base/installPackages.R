#!/usr/bin/env Rscript

install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "readr", "tximport", 
                        "optparse", "dplyr", "rtracklayer", "plyr", 
                        "survcomp", "VennDiagram", "knitr", "zoo",
                        "devtools", "plotly", "rmarkdown", "DT",
                        "heatmaply", "shiny"),ask = FALSE)

devtools::install_github("rstudio/d3heatmap", dependencies = TRUE, upgrade = "always")
remotes::install_github("alaimos/metaseqR", dependencies = TRUE, upgrade = "always")
