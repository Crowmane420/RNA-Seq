#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  defaults <- list(
    libdir = ".Rlibs",
    orgdb_package = "org.At.tair.db"
  )

  if (length(args) == 0) {
    return(defaults)
  }

  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) {
      stop("Arguments must use the form --name=value")
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=(.*)$", "\\1", arg)
    if (!key %in% names(defaults)) {
      stop(sprintf("Unknown argument: %s", key))
    }
    defaults[[key]] <- value
  }

  defaults
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(args$libdir, .libPaths()))

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "scales",
  "stringr",
  "tibble",
  "tidyr"
)

bioc_packages <- unique(c(
  "AnnotationDbi",
  "apeglm",
  "clusterProfiler",
  "DESeq2",
  "GENIE3",
  "GO.db",
  "igraph",
  "KEGGREST",
  "WGCNA",
  args$orgdb_package
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Dependency installation complete.")
