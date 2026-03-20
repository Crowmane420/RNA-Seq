#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  defaults <- list(
    counts = "counts.tsv",
    metadata = "sample_metadata.csv",
    outdir = "rnaseq_analysis",
    libdir = ".Rlibs",
    gene_id_col = "gene_id",
    transcript_id_col = "transcript_id(s)",
    sample_id_col = "sample_id",
    factor1_col = "genotype",
    factor2_col = "light_condition",
    factor1_ref = "",
    factor2_ref = "",
    replicate_col = "replicate",
    group_col = "",
    orgdb_package = "org.At.tair.db",
    gene_id_keytype = "TAIR",
    target_regex = "FTSH|FtsH",
    target_label = "FTSH",
    kegg_organism = "ath",
    min_count = "10",
    min_samples = "2"
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
args$transcript_id_col <- if (args$transcript_id_col %in% c("", "NA", "NULL")) "" else args$transcript_id_col
args$replicate_col <- if (args$replicate_col %in% c("", "NA", "NULL")) "" else args$replicate_col
args$group_col <- if (args$group_col %in% c("", "NA", "NULL")) "" else args$group_col
args$kegg_organism <- if (args$kegg_organism %in% c("", "NA", "NULL")) "" else args$kegg_organism
args$min_count <- as.integer(args$min_count)
args$min_samples <- as.integer(args$min_samples)
dir.create(args$libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(args$libdir, .libPaths()))

required_packages <- unique(c(
  "AnnotationDbi",
  "clusterProfiler",
  "DESeq2",
  "dplyr",
  "GENIE3",
  "ggplot2",
  "ggrepel",
  "GO.db",
  "igraph",
  "KEGGREST",
  "patchwork",
  "scales",
  "stringr",
  "tibble",
  "tidyr",
  "WGCNA",
  args$orgdb_package
))

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(sprintf("Missing required packages: %s", paste(missing_packages, collapse = ", ")))
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(DESeq2)
  library(dplyr)
  library(GENIE3)
  library(ggplot2)
  library(ggrepel)
  library(GO.db)
  library(igraph)
  library(KEGGREST)
  library(patchwork)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(WGCNA)
  library(args$orgdb_package, character.only = TRUE)
})

OrgDb <- get(args$orgdb_package)

WGCNA::disableWGCNAThreads()
allowWGCNAThreads <- get("allowWGCNAThreads", envir = asNamespace("WGCNA"))
try(suppressWarnings(allowWGCNAThreads(nThreads = 2)), silent = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

make_dir <- function(...) {
  path <- file.path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

write_table <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, quote = TRUE)
}

save_plot_dual <- function(plot_obj, stem, width = 10, height = 8, dpi = 400) {
  png_path <- file.path(fig_png_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_pdf_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  ggsave(pdf_path, plot = plot_obj, width = width, height = height, units = "in", device = grDevices::pdf, bg = "white")
}

matrix_to_long <- function(mat, row_name = "row_id", col_name = "col_id", value_name = "value") {
  long_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(long_df) <- c(row_name, col_name, value_name)
  long_df
}

rescale_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled[scaled > 2.5] <- 2.5
  scaled[scaled < -2.5] <- -2.5
  scaled
}

parse_ratio <- function(x) {
  as.numeric(sub("/.*$", "", x)) / as.numeric(sub("^.*/", "", x))
}

sanitize_token <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

ordered_levels <- function(values, reference_level = "") {
  values <- unique(as.character(values))
  if (length(values) == 0) {
    stop("Cannot determine factor levels from an empty vector.")
  }
  if (identical(reference_level, "") || is.na(reference_level)) {
    reference_level <- values[[1]]
  }
  if (!reference_level %in% values) {
    stop(sprintf("Reference level '%s' was not found in the data.", reference_level))
  }
  c(reference_level, values[values != reference_level])
}

build_palette <- function(levels_vec, seed_colors) {
  if (length(levels_vec) <= length(seed_colors)) {
    stats::setNames(seed_colors[seq_along(levels_vec)], levels_vec)
  } else {
    stats::setNames(grDevices::colorRampPalette(seed_colors)(length(levels_vec)), levels_vec)
  }
}

make_group_id <- function(factor1_value, factor2_value) {
  paste(sanitize_token(factor1_value), sanitize_token(factor2_value), sep = "__")
}

contrast_label_map <- character()
pretty_contrast_label <- function(x) {
  labels <- unname(contrast_label_map[x])
  labels[is.na(labels)] <- x[is.na(labels)]
  labels
}

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#495057"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
)

regulator_cols <- c(
  `Transcription factor` = "#C8553D",
  Kinase = "#2A6F97",
  `Transcription factor/Kinase` = "#6A4C93",
  Target = "#2D6A4F",
  Hub = "#6C757D"
)

outdir <- args$outdir
fig_png_dir <- make_dir(outdir, "figures", "png")
fig_pdf_dir <- make_dir(outdir, "figures", "pdf")
tables_dir <- make_dir(outdir, "tables")
tables_kegg_dir <- make_dir(tables_dir, "kegg")
tables_wgcna_dir <- make_dir(tables_dir, "wgcna")
tables_grn_dir <- make_dir(tables_dir, "grn")
tables_reg_dir <- make_dir(tables_dir, "regulators")
logs_dir <- make_dir(outdir, "logs")

message("Reading counts, metadata, and existing DE results...")

count_df <- read.delim(args$counts, check.names = FALSE)
metadata <- read.csv(args$metadata)

if (!args$gene_id_col %in% colnames(count_df)) {
  stop(sprintf("Count matrix is missing the configured gene ID column: %s", args$gene_id_col))
}
if (!identical(args$transcript_id_col, "") && !args$transcript_id_col %in% colnames(count_df)) {
  stop(sprintf("Count matrix is missing the configured transcript column: %s", args$transcript_id_col))
}

required_metadata_cols <- c(args$sample_id_col, args$factor1_col, args$factor2_col)
if (!all(required_metadata_cols %in% colnames(metadata))) {
  stop(
    sprintf(
      "Metadata must contain the configured sample/factor columns: %s",
      paste(required_metadata_cols, collapse = ", ")
    )
  )
}
if (!identical(args$replicate_col, "") && !args$replicate_col %in% colnames(metadata)) {
  stop(sprintf("Metadata is missing the configured replicate column: %s", args$replicate_col))
}
if (!identical(args$group_col, "") && !args$group_col %in% colnames(metadata)) {
  stop(sprintf("Metadata is missing the configured group column: %s", args$group_col))
}

colnames(count_df)[colnames(count_df) == args$gene_id_col] <- "gene_id"
if (!identical(args$transcript_id_col, "")) {
  colnames(count_df)[colnames(count_df) == args$transcript_id_col] <- "transcript_id(s)"
} else {
  count_df[["transcript_id(s)"]] <- NA_character_
}

colnames(metadata)[colnames(metadata) == args$sample_id_col] <- "sample_id"
colnames(metadata)[colnames(metadata) == args$factor1_col] <- "genotype"
colnames(metadata)[colnames(metadata) == args$factor2_col] <- "light_condition"
if (!identical(args$replicate_col, "")) {
  colnames(metadata)[colnames(metadata) == args$replicate_col] <- "replicate"
} else {
  metadata$replicate <- seq_len(nrow(metadata))
}
if (!identical(args$group_col, "")) {
  colnames(metadata)[colnames(metadata) == args$group_col] <- "group"
}

sample_columns <- setdiff(colnames(count_df), c("gene_id", "transcript_id(s)"))
metadata <- metadata[match(sample_columns, metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id
factor1_levels <- ordered_levels(metadata$genotype, args$factor1_ref)
factor2_levels <- ordered_levels(metadata$light_condition, args$factor2_ref)
if (length(factor2_levels) != 2) {
  stop("This reusable extension currently supports exactly two levels for factor2.")
}
metadata$genotype <- factor(metadata$genotype, levels = factor1_levels)
metadata$light_condition <- factor(metadata$light_condition, levels = factor2_levels)
metadata$group <- make_group_id(metadata$genotype, metadata$light_condition)
expected_groups <- apply(
  expand.grid(factor1_levels, factor2_levels, stringsAsFactors = FALSE),
  1,
  function(x) make_group_id(x[[1]], x[[2]])
)
metadata$group <- factor(metadata$group, levels = expected_groups)
metadata$sample_label <- paste(metadata$group, metadata$replicate, sep = "_")

factor1_name <- args$factor1_col
factor2_name <- args$factor2_col
factor1_ref <- factor1_levels[[1]]
factor2_ref <- factor2_levels[[1]]
factor2_case <- factor2_levels[[2]]
genotype_cols <- build_palette(factor1_levels, c("#3A7CA5", "#D1495B", "#2A9D8F", "#F4A261", "#6D597A"))
target_stem <- sanitize_token(tolower(args$target_label))
regulator_cols[[args$target_label]] <- regulator_cols[["Target"]]
regulator_cols <- regulator_cols[names(regulator_cols) != "Target"]

count_matrix <- as.matrix(count_df[, sample_columns, drop = FALSE])
rownames(count_matrix) <- count_df$gene_id
storage.mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ genotype * light_condition
)
dds <- dds[rowSums(counts(dds) >= args$min_count) >= args$min_samples, ]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

annotation_path <- file.path(outdir, "tables", "annotation", "gene_annotation.csv")
if (file.exists(annotation_path)) {
  gene_annotation <- read.csv(annotation_path)
} else {
  annotation_raw <- AnnotationDbi::select(
    OrgDb,
    keys = rownames(dds),
    keytype = args$gene_id_keytype,
    columns = c(args$gene_id_keytype, "SYMBOL", "GENENAME")
  )
  gene_annotation <- annotation_raw %>%
    transmute(gene_id = .data[[args$gene_id_keytype]], symbol = SYMBOL, gene_name = GENENAME) %>%
    group_by(gene_id) %>%
    summarise(
      symbol = first(na.omit(symbol)) %||% gene_id[[1]],
      gene_name = first(na.omit(gene_name)) %||% "",
      .groups = "drop"
    )
}

gene_annotation <- gene_annotation %>%
  mutate(
    symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol),
    gene_name = ifelse(is.na(gene_name), "", gene_name)
  )

symbol_lookup <- setNames(gene_annotation$symbol, gene_annotation$gene_id)

message("Building KEGG pathway annotations...")

kegg_links_raw <- if (!identical(args$kegg_organism, "")) {
  tryCatch(KEGGREST::keggLink("pathway", args$kegg_organism), error = function(e) NULL)
} else {
  NULL
}
kegg_names_raw <- if (!identical(args$kegg_organism, "")) {
  tryCatch(KEGGREST::keggList("pathway", args$kegg_organism), error = function(e) NULL)
} else {
  NULL
}

if (is.null(kegg_links_raw) || length(kegg_links_raw) == 0 || is.null(kegg_names_raw) || length(kegg_names_raw) == 0) {
  kegg_term2gene <- tibble(pathway_id = character(), gene_id = character())
  kegg_term2name <- tibble(pathway_id = character(), Description = character())
} else {
  kegg_term2gene <- tibble(
    gene_id = sub(paste0("^", args$kegg_organism, ":"), "", names(kegg_links_raw)),
    pathway_id = sub("^path:", "", unname(kegg_links_raw))
  ) %>%
    filter(gene_id %in% rownames(dds)) %>%
    distinct(pathway_id, gene_id)

  kegg_term2name <- tibble(
    pathway_id = names(kegg_names_raw),
    Description = sub(" - .*?$", "", unname(kegg_names_raw))
  ) %>%
    distinct(pathway_id, .keep_all = TRUE)
}

ftsh_direct_kegg <- tibble()

results_files <- list.files(file.path(outdir, "tables", "de_full"), pattern = "_full_results\\.csv$", full.names = TRUE)
if (length(results_files) == 0) {
  stop("No DE result tables found in the existing output directory.")
}

results_list <- setNames(
  lapply(results_files, function(path) read.csv(path, check.names = FALSE)),
  sub("_full_results\\.csv$", "", basename(results_files))
)

contrast_label_map <- stats::setNames(
  vapply(results_list, function(df) unique(df$contrast_label)[[1]], character(1)),
  names(results_list)
)

primary_mutant_contrasts <- unlist(lapply(factor2_levels, function(factor2_level) {
  vapply(setdiff(factor1_levels, factor1_ref), function(factor1_level) {
    sanitize_token(paste("factor1", factor1_level, "vs", factor1_ref, "in", factor2_name, factor2_level, sep = "_"))
  }, character(1))
}), use.names = FALSE)
factor2_contrasts <- vapply(factor1_levels, function(factor1_level) {
  sanitize_token(paste("factor2", factor2_case, "vs", factor2_ref, "in", factor1_name, factor1_level, sep = "_"))
}, character(1))
interaction_contrasts <- vapply(setdiff(factor1_levels, factor1_ref), function(factor1_level) {
  sanitize_token(paste("interaction", factor1_level, "vs", factor1_ref, factor2_name, "response", sep = "_"))
}, character(1))
kegg_gsea_targets <- c(primary_mutant_contrasts, interaction_contrasts)

message("Annotating transcription factors and kinases...")

get_go_descendants <- function(go_id) {
  descendants <- GO.db::GOMFOFFSPRING[[go_id]]
  unique(c(go_id, descendants[!is.na(descendants)]))
}

go_map <- AnnotationDbi::select(
  OrgDb,
  keys = rownames(dds),
  keytype = args$gene_id_keytype,
  columns = c(args$gene_id_keytype, "GOALL", "ONTOLOGYALL")
)

go_map <- go_map %>%
  filter(!is.na(GOALL), ONTOLOGYALL == "MF") %>%
  distinct(.data[[args$gene_id_keytype]], GOALL)

tf_terms <- unique(c(get_go_descendants("GO:0003700"), get_go_descendants("GO:0140110")))
kinase_terms <- unique(c(get_go_descendants("GO:0004672"), get_go_descendants("GO:0016301")))

tf_ids <- unique(go_map[[args$gene_id_keytype]][go_map$GOALL %in% tf_terms])
kinase_ids <- unique(go_map[[args$gene_id_keytype]][go_map$GOALL %in% kinase_terms])

regulator_annotation <- bind_rows(
  tibble(gene_id = tf_ids, regulator_class = "Transcription factor"),
  tibble(gene_id = kinase_ids, regulator_class = "Kinase")
) %>%
  filter(gene_id %in% rownames(vst_mat)) %>%
  distinct() %>%
  group_by(gene_id) %>%
  summarise(regulator_class = paste(sort(unique(regulator_class)), collapse = "/"), .groups = "drop") %>%
  left_join(gene_annotation, by = "gene_id") %>%
  arrange(regulator_class, symbol, gene_id)

write_table(regulator_annotation, file.path(tables_reg_dir, "regulator_annotation.csv"))

message("Running KEGG enrichment...")

run_kegg_ora <- function(gene_ids, label, direction, context) {
  if (nrow(kegg_term2gene) == 0 || nrow(kegg_term2name) == 0) {
    return(NULL)
  }
  gene_ids <- unique(gene_ids[gene_ids %in% rownames(dds)])
  if (length(gene_ids) < 15) {
    return(NULL)
  }

  kegg_obj <- tryCatch(
    clusterProfiler::enricher(
      gene = gene_ids,
      universe = rownames(dds),
      TERM2GENE = kegg_term2gene,
      TERM2NAME = kegg_term2name,
      minGSSize = 10,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    ),
    error = function(e) NULL
  )

  if (is.null(kegg_obj)) {
    return(NULL)
  }

  kegg_df <- as.data.frame(kegg_obj)
  if (nrow(kegg_df) == 0) {
    return(NULL)
  }

  kegg_df %>%
    mutate(
      label = label,
      direction = direction,
      context = context,
      gene_ratio_numeric = parse_ratio(GeneRatio)
    ) %>%
    arrange(p.adjust, desc(Count))
}

run_kegg_gsea <- function(res_df, contrast_name) {
  if (nrow(kegg_term2gene) == 0 || nrow(kegg_term2name) == 0) {
    return(NULL)
  }
  gene_stats <- res_df %>%
    filter(!is.na(stat), is.finite(stat)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(desc(stat))

  ranked_list <- gene_stats$stat
  names(ranked_list) <- gene_stats$gene_id
  ranked_list <- sort(ranked_list, decreasing = TRUE)

  gsea_obj <- tryCatch(
    clusterProfiler::GSEA(
      geneList = ranked_list,
      TERM2GENE = kegg_term2gene,
      TERM2NAME = kegg_term2name,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      seed = TRUE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(gsea_obj)) {
    return(NULL)
  }

  gsea_df <- as.data.frame(gsea_obj)
  if (nrow(gsea_df) == 0) {
    return(NULL)
  }

  gsea_df %>%
    mutate(
      contrast = contrast_name,
      contrast_label = pretty_contrast_label(contrast_name)
    ) %>%
    arrange(p.adjust, desc(abs(NES)))
}

kegg_ora_results <- list()
for (name in names(results_list)) {
  res_df <- results_list[[name]]
  kegg_ora_results[[paste0(name, "_up")]] <- run_kegg_ora(
    res_df$gene_id[res_df$significance == "Up-regulated"],
    pretty_contrast_label(name),
    "Up-regulated",
    "contrast"
  )
  kegg_ora_results[[paste0(name, "_down")]] <- run_kegg_ora(
    res_df$gene_id[res_df$significance == "Down-regulated"],
    pretty_contrast_label(name),
    "Down-regulated",
    "contrast"
  )
}

kegg_ora_summary <- bind_rows(kegg_ora_results)
if (nrow(kegg_ora_summary) == 0) {
  kegg_ora_summary <- tibble()
}

kegg_gsea_summary <- bind_rows(lapply(kegg_gsea_targets, function(name) run_kegg_gsea(results_list[[name]], name)))
if (nrow(kegg_gsea_summary) == 0) {
  kegg_gsea_summary <- tibble()
}

write_table(kegg_ora_summary, file.path(tables_kegg_dir, "kegg_ora_summary.csv"))
write_table(kegg_gsea_summary, file.path(tables_kegg_dir, "kegg_gsea_summary.csv"))

message("Running WGCNA...")

ftsh_ids <- gene_annotation %>%
  filter(
    grepl(args$target_regex, symbol, ignore.case = TRUE) |
      grepl(args$target_regex, gene_name, ignore.case = TRUE)
  ) %>%
  pull(gene_id) %>%
  unique() %>%
  sort()

gene_variance <- apply(vst_mat, 1, var)
base_wgcna_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(4000, length(gene_variance)))]
wgcna_genes <- unique(c(base_wgcna_genes, ftsh_ids, regulator_annotation$gene_id))
wgcna_genes <- intersect(wgcna_genes, rownames(vst_mat))

datExpr <- t(vst_mat[wgcna_genes, , drop = FALSE])
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

metadata_wgcna <- metadata[rownames(datExpr), , drop = FALSE]
powers <- c(1:10, seq(12, 20, 2))
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
fit_df <- as_tibble(sft$fitIndices)
soft_power <- fit_df %>%
  filter(SFT.R.sq >= 0.8) %>%
  summarise(power = min(Power)) %>%
  pull(power)
if (length(soft_power) == 0 || is.infinite(soft_power) || is.na(soft_power)) {
  soft_power <- fit_df$Power[which.max(fit_df$SFT.R.sq)]
}

net <- WGCNA::blockwiseModules(
  datExpr,
  power = soft_power,
  maxBlockSize = ncol(datExpr),
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 0
)

module_colors <- WGCNA::labels2colors(net$colors)
module_eigengenes <- WGCNA::orderMEs(net$MEs)
me_numeric_labels <- suppressWarnings(as.numeric(sub("^ME", "", colnames(module_eigengenes))))
me_color_labels <- ifelse(is.na(me_numeric_labels), sub("^ME", "", colnames(module_eigengenes)), WGCNA::labels2colors(me_numeric_labels))
colnames(module_eigengenes) <- paste0("ME", me_color_labels)

trait_df <- tibble(row.names = rownames(metadata_wgcna))
for (factor1_level in setdiff(factor1_levels, factor1_ref)) {
  trait_df[[sanitize_token(factor1_level)]] <- as.numeric(metadata_wgcna$genotype == factor1_level)
}
trait_df[[sanitize_token(factor2_case)]] <- as.numeric(metadata_wgcna$light_condition == factor2_case)
for (factor1_level in setdiff(factor1_levels, factor1_ref)) {
  combo_name <- sanitize_token(paste(factor1_level, factor2_case, sep = "_"))
  trait_df[[combo_name]] <- as.numeric(
    metadata_wgcna$genotype == factor1_level & metadata_wgcna$light_condition == factor2_case
  )
}
rownames(trait_df) <- rownames(metadata_wgcna)

module_trait_cor <- cor(module_eigengenes, trait_df, use = "pairwise.complete.obs")
module_trait_p <- WGCNA::corPvalueStudent(module_trait_cor, nSamples = nrow(datExpr))

kme_df <- as.data.frame(WGCNA::signedKME(datExpr, module_eigengenes))
kme_df$gene_id <- colnames(datExpr)

wgcna_gene_modules <- tibble(
  gene_id = colnames(datExpr),
  module_label = as.character(net$colors),
  module = module_colors
) %>%
  left_join(gene_annotation, by = "gene_id") %>%
  left_join(kme_df, by = "gene_id") %>%
  rowwise() %>%
  mutate(
    module_membership = {
      col_name <- paste0("kME", module)
      row_data <- dplyr::pick(dplyr::everything())
      if (col_name %in% names(row_data)) row_data[[col_name]][[1]] else NA_real_
    }
  ) %>%
  ungroup()

ftsh_module_assignments <- wgcna_gene_modules %>%
  filter(gene_id %in% ftsh_ids) %>%
  arrange(module, symbol, gene_id)
ftsh_modules <- unique(ftsh_module_assignments$module[ftsh_module_assignments$module != "grey"])

ftsh_module_hubs <- wgcna_gene_modules %>%
  filter(module %in% ftsh_modules) %>%
  arrange(module, desc(abs(module_membership))) %>%
  group_by(module) %>%
  slice_head(n = 25) %>%
  ungroup()

write_table(fit_df, file.path(tables_wgcna_dir, "wgcna_soft_threshold_table.csv"))
write_table(
  matrix_to_long(module_trait_cor, "module", "trait", "correlation") %>%
    left_join(matrix_to_long(module_trait_p, "module", "trait", "p_value"), by = c("module", "trait")) %>%
    mutate(module = sub("^ME", "", module)),
  file.path(tables_wgcna_dir, "wgcna_module_trait_correlations.csv")
)
write_table(wgcna_gene_modules, file.path(tables_wgcna_dir, "wgcna_gene_modules.csv"))
write_table(ftsh_module_assignments, file.path(tables_wgcna_dir, paste0(target_stem, "_module_assignments.csv")))
write_table(ftsh_module_hubs, file.path(tables_wgcna_dir, paste0(target_stem, "_module_hubs.csv")))

message(sprintf("Running %s-module KEGG enrichment...", args$target_label))

ftsh_direct_kegg <- kegg_term2gene %>%
  filter(gene_id %in% ftsh_ids) %>%
  left_join(gene_annotation, by = "gene_id") %>%
  left_join(kegg_term2name, by = "pathway_id") %>%
  arrange(Description, symbol, gene_id)
write_table(ftsh_direct_kegg, file.path(tables_kegg_dir, paste0(target_stem, "_direct_kegg_membership.csv")))

ftsh_module_kegg <- bind_rows(lapply(ftsh_modules, function(module_name) {
  module_genes <- wgcna_gene_modules$gene_id[wgcna_gene_modules$module == module_name]
  run_kegg_ora(module_genes, paste("Module", module_name), "Module genes", "ftsh_module")
}))
if (nrow(ftsh_module_kegg) == 0) {
  ftsh_module_kegg <- tibble()
}
write_table(ftsh_module_kegg, file.path(tables_kegg_dir, paste0(target_stem, "_module_kegg.csv")))

message("Building regulator summaries and GRN...")

regulator_deg_summary <- bind_rows(lapply(names(results_list), function(name) {
  res_df <- results_list[[name]]
  bind_rows(
    tibble(
      contrast = name,
      contrast_label = pretty_contrast_label(name),
      regulator_class = "Transcription factor",
      direction = "Up-regulated",
      count = sum(res_df$gene_id %in% tf_ids & res_df$significance == "Up-regulated", na.rm = TRUE)
    ),
    tibble(
      contrast = name,
      contrast_label = pretty_contrast_label(name),
      regulator_class = "Transcription factor",
      direction = "Down-regulated",
      count = sum(res_df$gene_id %in% tf_ids & res_df$significance == "Down-regulated", na.rm = TRUE)
    ),
    tibble(
      contrast = name,
      contrast_label = pretty_contrast_label(name),
      regulator_class = "Kinase",
      direction = "Up-regulated",
      count = sum(res_df$gene_id %in% kinase_ids & res_df$significance == "Up-regulated", na.rm = TRUE)
    ),
    tibble(
      contrast = name,
      contrast_label = pretty_contrast_label(name),
      regulator_class = "Kinase",
      direction = "Down-regulated",
      count = sum(res_df$gene_id %in% kinase_ids & res_df$significance == "Down-regulated", na.rm = TRUE)
    )
  )
}))

write_table(regulator_deg_summary, file.path(tables_reg_dir, "regulator_deg_summary.csv"))

regulator_lfc_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), lapply(names(results_list), function(name) {
  results_list[[name]] %>%
    select(gene_id, log2FoldChange) %>%
    rename(!!name := log2FoldChange)
}))

regulator_rank <- bind_rows(lapply(primary_mutant_contrasts, function(name) {
  results_list[[name]] %>%
    select(gene_id, log2FoldChange)
})) %>%
  group_by(gene_id) %>%
  summarise(max_abs_log2FC = max(abs(log2FoldChange), na.rm = TRUE), .groups = "drop")

top_tf_ids <- regulator_rank %>%
  filter(gene_id %in% tf_ids) %>%
  arrange(desc(max_abs_log2FC)) %>%
  slice_head(n = 25) %>%
  pull(gene_id)

top_kinase_ids <- regulator_rank %>%
  filter(gene_id %in% kinase_ids) %>%
  arrange(desc(max_abs_log2FC)) %>%
  slice_head(n = 25) %>%
  pull(gene_id)

write_table(regulator_rank %>% filter(gene_id %in% top_tf_ids), file.path(tables_reg_dir, "top_transcription_factors.csv"))
write_table(regulator_rank %>% filter(gene_id %in% top_kinase_ids), file.path(tables_reg_dir, "top_kinases.csv"))

ftsh_module_regulators <- wgcna_gene_modules %>%
  filter(module %in% ftsh_modules, gene_id %in% regulator_annotation$gene_id) %>%
  arrange(desc(abs(module_membership)))

candidate_regulators <- unique(c(
  ftsh_module_regulators$gene_id,
  regulator_rank %>% filter(gene_id %in% regulator_annotation$gene_id) %>% arrange(desc(max_abs_log2FC)) %>% slice_head(n = 60) %>% pull(gene_id)
))

candidate_targets <- unique(c(
  ftsh_ids,
  ftsh_module_hubs %>% filter(!gene_id %in% candidate_regulators) %>% pull(gene_id)
))

network_genes <- unique(c(candidate_regulators, candidate_targets))
network_genes <- intersect(network_genes, rownames(vst_mat))
expr_network <- vst_mat[network_genes, , drop = FALSE]
expr_network <- expr_network[apply(expr_network, 1, var) > 0, , drop = FALSE]
network_regulators <- intersect(candidate_regulators, rownames(expr_network))
network_targets <- intersect(candidate_targets, rownames(expr_network))

if (length(network_regulators) >= 2 && nrow(expr_network) >= 3 && length(network_targets) >= 1) {
  weight_matrix <- GENIE3::GENIE3(expr_network, regulators = network_regulators, nTrees = 500, nCores = 1)
  link_list <- GENIE3::getLinkList(weight_matrix)
  link_list <- as_tibble(link_list) %>%
    mutate(
      regulatoryGene = as.character(regulatoryGene),
      targetGene = as.character(targetGene)
    ) %>%
    filter(regulatoryGene != targetGene, targetGene %in% network_targets) %>%
    arrange(desc(weight))
} else {
  link_list <- tibble(regulatoryGene = character(), targetGene = character(), weight = numeric())
}

if (nrow(link_list) > 120) {
  cutoff <- quantile(link_list$weight, probs = 0.98, na.rm = TRUE)
  filtered_links <- link_list %>% filter(weight >= cutoff)
  if (nrow(filtered_links) < 80) {
    filtered_links <- link_list %>% slice_head(n = 80)
  } else if (nrow(filtered_links) > 150) {
    filtered_links <- filtered_links %>% slice_head(n = 150)
  }
  link_list <- filtered_links
}

mutation_scores <- bind_rows(lapply(primary_mutant_contrasts, function(name) {
  results_list[[name]] %>% select(gene_id, log2FoldChange, padj)
})) %>%
  group_by(gene_id) %>%
  summarise(
    mutation_score = max(abs(log2FoldChange), na.rm = TRUE),
    min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(min_padj = ifelse(is.infinite(min_padj), NA_real_, min_padj))

grn_nodes <- tibble(gene_id = unique(c(link_list$regulatoryGene, link_list$targetGene))) %>%
  left_join(gene_annotation, by = "gene_id") %>%
  left_join(regulator_annotation %>% select(gene_id, regulator_class), by = "gene_id") %>%
  left_join(mutation_scores, by = "gene_id") %>%
  mutate(
    node_class = case_when(
      gene_id %in% ftsh_ids ~ args$target_label,
      !is.na(regulator_class) ~ regulator_class,
      TRUE ~ "Hub"
    )
  )

if (nrow(link_list) > 0 && nrow(grn_nodes) > 0) {
  graph_obj <- igraph::graph_from_data_frame(
    d = link_list %>% rename(from = regulatoryGene, to = targetGene),
    directed = TRUE,
    vertices = grn_nodes %>% rename(name = gene_id)
  )

  coords <- igraph::layout_with_fr(graph_obj, weights = E(graph_obj)$weight)
  coord_df <- tibble(
    gene_id = V(graph_obj)$name,
    x = coords[, 1],
    y = coords[, 2],
    degree = igraph::degree(graph_obj, mode = "all")
  )
} else {
  coord_df <- tibble(gene_id = character(), x = numeric(), y = numeric(), degree = numeric())
}

label_degree_cutoff <- if (nrow(coord_df) > 0) quantile(coord_df$degree, 0.85, na.rm = TRUE) else Inf

grn_nodes_plot <- grn_nodes %>%
  left_join(coord_df, by = "gene_id") %>%
  mutate(
    label = case_when(
      gene_id %in% ftsh_ids ~ symbol,
      degree >= label_degree_cutoff ~ symbol,
      TRUE ~ NA_character_
    ),
    point_size = ifelse(is.na(degree), 2.5, pmax(2.5, rescale(degree, to = c(2.5, 8))))
  )

grn_edges_plot <- link_list %>%
  rename(from = regulatoryGene, to = targetGene) %>%
  left_join(coord_df %>% rename(from = gene_id, x_from = x, y_from = y), by = "from") %>%
  left_join(coord_df %>% rename(to = gene_id, x_to = x, y_to = y), by = "to")

write_table(link_list, file.path(tables_grn_dir, paste0(target_stem, "_grn_edges.csv")))
write_table(grn_nodes_plot, file.path(tables_grn_dir, paste0(target_stem, "_grn_nodes.csv")))

message("Building figures...")

kegg_plot <- NULL
if (nrow(kegg_ora_summary) > 0) {
  kegg_plot_df <- kegg_ora_summary %>%
    filter(label %in% pretty_contrast_label(c(primary_mutant_contrasts, factor2_contrasts))) %>%
    group_by(label, direction) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    mutate(Description = stringr::str_wrap(Description, width = 36))

  kegg_plot <- ggplot(kegg_plot_df, aes(x = gene_ratio_numeric, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D1495B", high = "#2F6690", trans = "reverse") +
    facet_grid(direction ~ label, scales = "free_y", space = "free_y") +
    labs(
      title = "KEGG pathway enrichment across contrasts",
      subtitle = "Top over-represented KEGG terms among DEGs in each direction",
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Count"
    ) +
    theme(axis.text.y = element_text(size = 8), strip.text.x = element_text(size = 10))
  save_plot_dual(kegg_plot, "kegg_enrichment", width = 18, height = 14)
}

ftsh_module_kegg_plot <- NULL
if (nrow(ftsh_module_kegg) > 0) {
  ftsh_module_kegg_df <- ftsh_module_kegg %>%
    group_by(label) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(Description = stringr::str_wrap(Description, width = 36))

  ftsh_module_kegg_plot <- ggplot(ftsh_module_kegg_df, aes(x = gene_ratio_numeric, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D1495B", high = "#2F6690", trans = "reverse") +
    facet_wrap(~ label, scales = "free_y") +
    labs(
      title = sprintf("KEGG context of %s-associated WGCNA modules", args$target_label),
      subtitle = sprintf("Pathways enriched among genes in modules containing %s-associated loci", args$target_label),
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Count"
    ) +
    theme(axis.text.y = element_text(size = 8))
  save_plot_dual(ftsh_module_kegg_plot, paste0(target_stem, "_module_kegg"), width = 16, height = 10)
}

soft_threshold_plot <- ggplot(fit_df, aes(Power, SFT.R.sq)) +
  geom_line(color = "#2F6690", linewidth = 0.7) +
  geom_point(color = "#2F6690", size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#D1495B") +
  geom_vline(xintercept = soft_power, linetype = "dotted", color = "#2D6A4F") +
  labs(
    title = "WGCNA soft-threshold selection",
    subtitle = sprintf("Chosen power = %s", soft_power),
    x = "Soft-thresholding power",
    y = "Scale-free topology fit (signed R^2)"
  )

module_trait_df <- matrix_to_long(module_trait_cor, "module", "trait", "correlation") %>%
  left_join(matrix_to_long(module_trait_p, "module", "trait", "p_value"), by = c("module", "trait")) %>%
  mutate(
    module = sub("^ME", "", module),
    module = factor(module, levels = rev(sub("^ME", "", colnames(module_eigengenes)))),
    label = sprintf("%.2f\np=%s", correlation, signif(p_value, 2))
  )

module_trait_plot <- ggplot(module_trait_df, aes(trait, module, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3.1) +
  scale_fill_gradient2(low = "#2F6690", mid = "white", high = "#D1495B", midpoint = 0) +
  labs(
    title = "Module-trait correlations",
    subtitle = sprintf("Exploratory WGCNA run on highly variable genes plus %s/regulator loci", args$target_label),
    x = NULL,
    y = NULL,
    fill = "Correlation"
  )

ftsh_module_eigengene_plot <- NULL
if (length(ftsh_modules) > 0) {
  ftsh_module_me_names <- paste0("ME", ftsh_modules)
  ftsh_module_me_names <- ftsh_module_me_names[ftsh_module_me_names %in% colnames(module_eigengenes)]
  if (length(ftsh_module_me_names) > 0) {
    eigengene_long <- as.data.frame(module_eigengenes[, ftsh_module_me_names, drop = FALSE]) %>%
      rownames_to_column("sample_id") %>%
      left_join(metadata_wgcna %>% as.data.frame() %>% tibble::as_tibble(), by = "sample_id") %>%
      pivot_longer(cols = all_of(ftsh_module_me_names), names_to = "module", values_to = "eigengene") %>%
      mutate(module = sub("^ME", "", module))

    ftsh_module_eigengene_plot <- ggplot(eigengene_long, aes(group, eigengene, color = genotype, shape = light_condition)) +
      geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.4, alpha = 0.9) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 2.8, color = "black") +
      scale_color_manual(values = genotype_cols) +
      facet_wrap(~ module, scales = "free_y") +
      labs(
      title = sprintf("%s-module eigengene behavior across groups", args$target_label),
      subtitle = sprintf("Modules shown here contain at least one %s-associated gene", args$target_label),
        x = NULL,
        y = "Module eigengene",
        color = "Genotype",
        shape = "Light"
      ) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save_plot_dual(ftsh_module_eigengene_plot, paste0("wgcna_", target_stem, "_modules"), width = 16, height = 10)
  }
}

wgcna_composite <- soft_threshold_plot / module_trait_plot +
  plot_annotation(
    title = "WGCNA overview",
    subtitle = sprintf("Coexpression modules were used to connect %s and %s effects with %s-associated transcriptional programs.", factor1_name, factor2_name, args$target_label),
    tag_levels = "A"
  )
save_plot_dual(wgcna_composite, "wgcna_overview", width = 14, height = 14)

reg_counts_plot <- regulator_deg_summary %>%
  mutate(
    signed_count = ifelse(direction == "Down-regulated", -count, count),
    contrast_label = factor(contrast_label, levels = pretty_contrast_label(names(results_list)))
  ) %>%
  ggplot(aes(contrast_label, signed_count, fill = direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  facet_wrap(~ regulator_class, scales = "free_y") +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c(`Up-regulated` = "#D1495B", `Down-regulated` = "#2F6690")) +
  labs(
    title = "Differentially expressed regulators by contrast",
    x = NULL,
    y = "Count",
    fill = NULL
  )

make_regulator_heatmap <- function(gene_ids, title) {
  if (length(gene_ids) == 0) {
    return(NULL)
  }
  contrast_order <- c(primary_mutant_contrasts, factor2_contrasts, interaction_contrasts)
  heat_df <- bind_rows(lapply(contrast_order, function(name) {
    results_list[[name]] %>%
      filter(gene_id %in% gene_ids) %>%
      transmute(gene_id, contrast = pretty_contrast_label(name), log2FoldChange)
  })) %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol))

  gene_levels <- regulator_rank %>%
    filter(gene_id %in% gene_ids) %>%
    arrange(desc(max_abs_log2FC)) %>%
    pull(gene_id)

  heat_df %>%
    mutate(
      symbol = factor(symbol, levels = rev(unique(gene_annotation$symbol[match(gene_levels, gene_annotation$gene_id)]))),
      contrast = factor(contrast, levels = pretty_contrast_label(contrast_order))
    ) %>%
    ggplot(aes(contrast, symbol, fill = log2FoldChange)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2F6690", mid = "white", high = "#D1495B", midpoint = 0) +
    labs(title = title, x = NULL, y = NULL, fill = "log2FC") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 7))
}

tf_heatmap <- make_regulator_heatmap(top_tf_ids, "Top transcription factors")
kinase_heatmap <- make_regulator_heatmap(top_kinase_ids, "Top kinases")
if (!is.null(tf_heatmap)) {
  save_plot_dual(tf_heatmap, "transcription_factor_heatmap", width = 12, height = 10)
}
if (!is.null(kinase_heatmap)) {
  save_plot_dual(kinase_heatmap, "kinase_heatmap", width = 12, height = 10)
}
regulator_plot_list <- c(list(reg_counts_plot), Filter(Negate(is.null), list(tf_heatmap, kinase_heatmap)))
regulator_composite <- wrap_plots(
  regulator_plot_list,
  ncol = 1,
  heights = c(1, rep(1.2, length(regulator_plot_list) - 1))
) +
  plot_annotation(
    title = "Regulator-focused view of the experiment",
    subtitle = "Separate summaries for transcription factors and kinases across mutant and light contrasts.",
    tag_levels = "A"
  )
save_plot_dual(regulator_composite, "regulator_overview", width = 16, height = 20)

grn_plot <- if (nrow(grn_edges_plot) > 0 && nrow(grn_nodes_plot) > 0) {
  ggplot() +
    geom_segment(
      data = grn_edges_plot,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = weight),
      color = "#8D99AE",
      linewidth = 0.5,
      arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
    ) +
    geom_point(
      data = grn_nodes_plot,
      aes(x = x, y = y, color = node_class, size = point_size),
      alpha = 0.95
    ) +
    ggrepel::geom_text_repel(
      data = subset(grn_nodes_plot, !is.na(label)),
      aes(x = x, y = y, label = label, color = node_class),
      size = 3.2,
      box.padding = 0.3,
      max.overlaps = 40,
      show.legend = FALSE
    ) +
    scale_color_manual(values = regulator_cols) +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.2, 0.9)) +
    labs(
      title = sprintf("%s-centered exploratory regulatory network", args$target_label),
      subtitle = sprintf("Edges were inferred with GENIE3 using candidate transcription factors, kinases, %s genes, and %s-module hubs.", args$target_label, args$target_label),
      x = NULL,
      y = NULL,
      color = NULL,
      alpha = "GENIE3 weight"
    ) +
    theme_void() +
    theme(legend.position = "bottom")
} else {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No GRN edges passed the filtering threshold for this focused network.",
      size = 5,
      color = "#495057"
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(
      title = sprintf("%s-centered exploratory regulatory network", args$target_label),
      subtitle = "GENIE3 did not yield a stable edge set after filtering in this focused network run.",
      x = NULL,
      y = NULL
    ) +
    theme_void()
}
save_plot_dual(grn_plot, paste0(target_stem, "_grn_network"), width = 14, height = 11)

message("Writing extension report...")

top_kegg_lines <- if (nrow(kegg_ora_summary) > 0) {
  kegg_ora_summary %>%
    group_by(label, direction) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(line = sprintf("- %s (%s): %s [padj = %.3g]", label, direction, Description, p.adjust)) %>%
    pull(line)
} else {
  "- No KEGG ORA term met the reporting threshold."
}

top_ftsh_module_lines <- if (nrow(ftsh_module_kegg) > 0) {
  ftsh_module_kegg %>%
    group_by(label) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(line = sprintf("- %s: %s [padj = %.3g]", label, Description, p.adjust)) %>%
    pull(line)
} else {
  sprintf("- No %s-associated module reached the KEGG reporting cutoff.", args$target_label)
}

ftsh_direct_lines <- if (nrow(ftsh_direct_kegg) > 0) {
  ftsh_direct_kegg %>%
    mutate(symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol)) %>%
    group_by(Description) %>%
    summarise(genes = paste(sort(unique(symbol)), collapse = ", "), .groups = "drop") %>%
    transmute(line = sprintf("- %s: %s", Description, genes)) %>%
    slice_head(n = 8) %>%
    pull(line)
} else {
  sprintf("- No %s-family gene had a direct KEGG pathway annotation.", args$target_label)
}

top_grn_lines <- link_list %>%
  left_join(gene_annotation %>% select(gene_id, symbol), by = c("regulatoryGene" = "gene_id")) %>%
  rename(regulator_symbol = symbol) %>%
  left_join(gene_annotation %>% select(gene_id, symbol), by = c("targetGene" = "gene_id")) %>%
  rename(target_symbol = symbol) %>%
  mutate(
    regulator_symbol = ifelse(is.na(regulator_symbol), regulatoryGene, regulator_symbol),
    target_symbol = ifelse(is.na(target_symbol), targetGene, target_symbol),
    line = sprintf("- %s -> %s [weight = %.4f]", regulator_symbol, target_symbol, weight)
  ) %>%
  slice_head(n = 10) %>%
  pull(line)

if (length(top_grn_lines) == 0) {
  top_grn_lines <- "- No GRN edges passed the reporting threshold."
}

kegg_md_line <- if (nrow(kegg_ora_summary) > 0) "![KEGG enrichment](figures/png/kegg_enrichment.png)" else NULL
ftsh_module_md_line <- if (nrow(ftsh_module_kegg) > 0) sprintf("![%s-module KEGG](figures/png/%s_module_kegg.png)", args$target_label, target_stem) else NULL
ftsh_module_eigengene_md_line <- if (!is.null(ftsh_module_eigengene_plot)) sprintf("![%s-module eigengenes](figures/png/wgcna_%s_modules.png)", args$target_label, target_stem) else NULL

kegg_html_line <- if (nrow(kegg_ora_summary) > 0) "<img src='figures/png/kegg_enrichment.png' alt='KEGG enrichment' />" else NULL
ftsh_module_html_line <- if (nrow(ftsh_module_kegg) > 0) sprintf("<img src='figures/png/%s_module_kegg.png' alt='%s module KEGG' />", target_stem, args$target_label) else NULL
ftsh_module_eigengene_html_line <- if (!is.null(ftsh_module_eigengene_plot)) sprintf("<img src='figures/png/wgcna_%s_modules.png' alt='%s module eigengenes' />", target_stem, args$target_label) else NULL

report_md <- c(
  "# Advanced pathway and network extension",
  "",
  "## Scope note",
  "",
  "- This extension adds KEGG, WGCNA, exploratory GRN inference, and separate regulator views to the existing RNA-seq analysis bundle.",
  "- WGCNA and GRN outputs should be interpreted as exploratory because the design includes 15 total samples and 2-3 replicates per group.",
  "",
  "## KEGG pathway enrichment",
  "",
  top_kegg_lines,
  "",
  kegg_md_line,
  "",
  sprintf("## %s-related pathway context", args$target_label),
  "",
  ftsh_direct_lines,
  "",
  top_ftsh_module_lines,
  "",
  ftsh_module_md_line,
  "",
  "## WGCNA",
  "",
  sprintf("- Soft-thresholding power selected for WGCNA: %s", soft_power),
  sprintf("- %s-associated modules detected: %s", args$target_label, if (length(ftsh_modules) > 0) paste(ftsh_modules, collapse = ", ") else "none outside grey"),
  "",
  "![WGCNA overview](figures/png/wgcna_overview.png)",
  "",
  ftsh_module_eigengene_md_line,
  "",
  "## Regulator-focused view",
  "",
  "- Separate summaries were generated for transcription factors and kinases using GO molecular-function annotations.",
  "",
  "![Regulator overview](figures/png/regulator_overview.png)",
  "",
  sprintf("## Exploratory %s-centered GRN", args$target_label),
  "",
  top_grn_lines,
  "",
  sprintf("![%s GRN](figures/png/%s_grn_network.png)", args$target_label, target_stem),
  ""
)

writeLines(report_md, file.path(outdir, "advanced_regulatory_report.md"))

report_html <- c(
  "<!DOCTYPE html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8' />",
  "<title>Advanced pathway and network extension</title>",
  "<style>",
  "body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; line-height: 1.6; color: #1f2933; padding: 0 20px; }",
  "h1, h2 { color: #102A43; }",
  "img { max-width: 100%; border: 1px solid #d9e2ec; margin: 16px 0 28px; }",
  "code { background: #f0f4f8; padding: 2px 4px; border-radius: 4px; }",
  "</style>",
  "</head>",
  "<body>",
  "<h1>Advanced pathway and network extension</h1>",
  "<h2>Scope note</h2>",
  "<ul>",
  "<li>This extension adds KEGG, WGCNA, exploratory GRN inference, and separate regulator views to the existing RNA-seq analysis bundle.</li>",
  "<li>WGCNA and GRN outputs should be interpreted as exploratory because the design includes 15 total samples and 2-3 replicates per group.</li>",
  "</ul>",
  "<h2>KEGG pathway enrichment</h2>",
  "<ul>",
  paste0("<li>", sub("^- ", "", top_kegg_lines), "</li>", collapse = ""),
  "</ul>",
  kegg_html_line,
  sprintf("<h2>%s-related pathway context</h2>", args$target_label),
  "<ul>",
  paste0("<li>", sub("^- ", "", ftsh_direct_lines), "</li>", collapse = ""),
  paste0("<li>", sub("^- ", "", top_ftsh_module_lines), "</li>", collapse = ""),
  "</ul>",
  ftsh_module_html_line,
  "<h2>WGCNA</h2>",
  sprintf("<p>Soft-thresholding power selected for WGCNA: <code>%s</code></p>", soft_power),
  sprintf("<p>%s-associated modules detected: <code>%s</code></p>", args$target_label, if (length(ftsh_modules) > 0) paste(ftsh_modules, collapse = ", ") else "none outside grey"),
  "<img src='figures/png/wgcna_overview.png' alt='WGCNA overview' />",
  ftsh_module_eigengene_html_line,
  "<h2>Regulator-focused view</h2>",
  "<p>Separate summaries were generated for transcription factors and kinases using GO molecular-function annotations.</p>",
  "<img src='figures/png/regulator_overview.png' alt='Regulator overview' />",
  sprintf("<h2>Exploratory %s-centered GRN</h2>", args$target_label),
  "<ul>",
  paste0("<li>", sub("^- ", "", top_grn_lines), "</li>", collapse = ""),
  "</ul>",
  sprintf("<img src='figures/png/%s_grn_network.png' alt='%s GRN' />", target_stem, args$target_label),
  "</body>",
  "</html>"
)

writeLines(report_html, file.path(outdir, "advanced_regulatory_report.html"))
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "advanced_extension_session_info.txt"))

message("Advanced pathway and network extension completed successfully.")
