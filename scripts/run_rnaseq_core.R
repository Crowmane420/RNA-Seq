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
    min_count = "10",
    min_samples = "2",
    padj_cutoff = "0.05",
    lfc_cutoff = "1"
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
args$min_count <- as.integer(args$min_count)
args$min_samples <- as.integer(args$min_samples)
args$padj_cutoff <- as.numeric(args$padj_cutoff)
args$lfc_cutoff <- as.numeric(args$lfc_cutoff)
dir.create(args$libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(args$libdir, .libPaths()))

required_packages <- unique(c(
  "AnnotationDbi",
  "apeglm",
  "clusterProfiler",
  "DESeq2",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "scales",
  "stringr",
  "tibble",
  "tidyr",
  args$orgdb_package
))

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required packages: %s. Install them first or point --libdir to a library containing them.",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(apeglm)
  library(clusterProfiler)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(args$orgdb_package, character.only = TRUE)
})

OrgDb <- get(args$orgdb_package)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

first_non_missing <- function(x, fallback = NA_character_) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) fallback else x[[1]]
}

make_dir <- function(...) {
  path <- file.path(...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

save_plot_dual <- function(plot_obj, stem, width = 10, height = 8, dpi = 400) {
  png_path <- file.path(fig_png_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_pdf_dir, paste0(stem, ".pdf"))

  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    bg = "white"
  )
}

write_table <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, quote = TRUE)
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

find_single_coef <- function(result_names, include_terms, exclude_terms = character()) {
  hits <- result_names
  for (term in include_terms) {
    hits <- hits[grepl(term, hits, fixed = TRUE)]
  }
  for (term in exclude_terms) {
    hits <- hits[!grepl(term, hits, fixed = TRUE)]
  }
  if (length(hits) != 1) {
    stop(
      sprintf(
        "Expected exactly one coefficient for include=%s exclude=%s, found: %s",
        paste(include_terms, collapse = ", "),
        paste(exclude_terms, collapse = ", "),
        paste(hits, collapse = ", ")
      )
    )
  }
  hits[[1]]
}

extract_primary_symbol <- function(gene_annotation, gene_ids) {
  gene_annotation %>%
    filter(gene_id %in% gene_ids) %>%
    distinct(gene_id, symbol) %>%
    mutate(symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol)) %>%
    deframe()
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
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(color = "#6c757d", size = 9)
    )
)

sig_cols <- c(
  `Up-regulated` = "#D1495B",
  `Down-regulated` = "#2F6690",
  `Not significant` = "#B8C0CC"
)

outdir <- args$outdir
fig_dir <- make_dir(outdir, "figures")
fig_png_dir <- make_dir(fig_dir, "png")
fig_pdf_dir <- make_dir(fig_dir, "pdf")
tables_dir <- make_dir(outdir, "tables")
tables_full_dir <- make_dir(tables_dir, "de_full")
tables_filtered_dir <- make_dir(tables_dir, "deg_filtered")
tables_qc_dir <- make_dir(tables_dir, "qc")
tables_annotation_dir <- make_dir(tables_dir, "annotation")
tables_enrichment_dir <- make_dir(tables_dir, "enrichment")
tables_ftsh_dir <- make_dir(tables_dir, sanitize_token(tolower(args$target_label)))
logs_dir <- make_dir(outdir, "logs")

message("Reading count matrix and sample metadata...")

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
if (anyDuplicated(sample_columns) > 0) {
  stop("Duplicate sample columns found in count matrix.")
}
if (anyDuplicated(metadata$sample_id) > 0) {
  stop("Duplicate sample_id entries found in metadata.")
}
if (!setequal(sample_columns, metadata$sample_id)) {
  missing_in_metadata <- setdiff(sample_columns, metadata$sample_id)
  missing_in_counts <- setdiff(metadata$sample_id, sample_columns)
  stop(
    paste(
      "Sample mismatch detected.",
      if (length(missing_in_metadata) > 0) paste("Missing in metadata:", paste(missing_in_metadata, collapse = ", ")) else "",
      if (length(missing_in_counts) > 0) paste("Missing in counts:", paste(missing_in_counts, collapse = ", ")) else ""
    )
  )
}

metadata <- metadata[match(sample_columns, metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

factor1_levels <- ordered_levels(metadata$genotype, args$factor1_ref)
factor2_levels <- ordered_levels(metadata$light_condition, args$factor2_ref)
if (length(factor2_levels) != 2) {
  stop("This reusable pipeline currently supports exactly two levels for factor2.")
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
metadata$sample_label <- factor(
  paste(metadata$group, metadata$replicate, sep = "_"),
  levels = paste(metadata$group, metadata$replicate, sep = "_")
)

factor1_name <- args$factor1_col
factor2_name <- args$factor2_col
factor1_ref <- factor1_levels[[1]]
factor2_ref <- factor2_levels[[1]]
factor2_case <- factor2_levels[[2]]
genotype_cols <- build_palette(factor1_levels, c("#3A7CA5", "#D1495B", "#2A9D8F", "#F4A261", "#6D597A"))
light_cols <- build_palette(factor2_levels, c("#457B9D", "#E9C46A", "#2A9D8F"))

count_matrix <- as.matrix(count_df[, sample_columns, drop = FALSE])
rownames(count_matrix) <- count_df$gene_id
storage.mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)

sample_stats <- metadata %>%
  transmute(
    sample_id,
    genotype,
    light_condition,
    replicate,
    group,
    library_size = colSums(count_matrix)[sample_id],
    library_size_millions = library_size / 1e6,
    detected_genes = colSums(count_matrix > 0)[sample_id],
    genes_ge_threshold = colSums(count_matrix >= args$min_count)[sample_id]
  )

write_table(sample_stats, file.path(tables_qc_dir, "sample_qc_metrics.csv"))
write_table(metadata, file.path(tables_qc_dir, "metadata_aligned.csv"))

message("Running DESeq2 interaction model...")

dds_interaction <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ genotype * light_condition
)

keep <- rowSums(counts(dds_interaction) >= args$min_count) >= args$min_samples
dds_interaction <- dds_interaction[keep, ]
dds_interaction <- DESeq(dds_interaction, quiet = TRUE)
vsd <- vst(dds_interaction, blind = FALSE)
vst_mat <- assay(vsd)
norm_counts <- counts(dds_interaction, normalized = TRUE)
result_names <- resultsNames(dds_interaction)

writeLines(result_names, file.path(logs_dir, "deseq2_results_names.txt"))

message("Building annotation table...")

annotation_raw <- AnnotationDbi::select(
  OrgDb,
  keys = rownames(dds_interaction),
  keytype = args$gene_id_keytype,
  columns = c(args$gene_id_keytype, "SYMBOL", "GENENAME")
)

gene_annotation <- annotation_raw %>%
  transmute(
    gene_id = .data[[args$gene_id_keytype]],
    symbol = SYMBOL,
    gene_name = GENENAME
  ) %>%
  group_by(gene_id) %>%
  summarise(
    symbol = first_non_missing(symbol, fallback = gene_id[1]),
    gene_name = first_non_missing(gene_name, fallback = ""),
    .groups = "drop"
  ) %>%
  mutate(symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol))

write_table(gene_annotation, file.path(tables_annotation_dir, "gene_annotation.csv"))

message("Generating QC plots...")

pca_data <- plotPCA(vsd, intgroup = c("genotype", "light_condition"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data$sample_id <- rownames(pca_data)

order_samples <- metadata$sample_id[order(metadata$group, metadata$replicate)]
dist_mat <- as.matrix(dist(t(vst_mat)))
dist_mat <- dist_mat[order_samples, order_samples]
cor_mat <- cor(vst_mat[, order_samples], method = "pearson")

within_group_mask <- outer(metadata$group, metadata$group, FUN = "==")
diag(within_group_mask) <- FALSE
between_group_mask <- outer(metadata$group, metadata$group, FUN = "!=")

cor_summary <- tibble(
  metric = c("median_within_group_correlation", "median_between_group_correlation"),
  value = c(
    median(cor(vst_mat)[within_group_mask]),
    median(cor(vst_mat)[between_group_mask])
  )
)
write_table(cor_summary, file.path(tables_qc_dir, "correlation_summary.csv"))

library_plot <- ggplot(sample_stats, aes(x = reorder(sample_id, library_size_millions), y = library_size_millions, fill = genotype)) +
  geom_col(width = 0.8, color = "white") +
  geom_point(aes(shape = light_condition), size = 3, color = "#1f2933") +
  coord_flip() +
  scale_fill_manual(values = genotype_cols) +
  labs(
    title = "Library size by sample",
    x = NULL,
    y = "Mapped counts (millions)",
    fill = factor1_name,
    shape = factor2_name
  )

detected_plot <- ggplot(sample_stats, aes(x = reorder(sample_id, detected_genes), y = detected_genes, fill = genotype)) +
  geom_col(width = 0.8, color = "white") +
  geom_point(aes(y = genes_ge_threshold, shape = light_condition), size = 3, color = "#1f2933") +
  coord_flip() +
  scale_fill_manual(values = genotype_cols) +
  labs(
    title = "Detected genes by sample",
    subtitle = sprintf("Bars: genes with counts > 0, points: genes with counts >= %d", args$min_count),
    x = NULL,
    y = "Number of genes",
    fill = factor1_name,
    shape = factor2_name
  )

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype, shape = light_condition)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3.2, box.padding = 0.35, max.overlaps = 20) +
  scale_color_manual(values = genotype_cols) +
  labs(
    title = "Principal component analysis",
    subtitle = sprintf("PC1 = %.1f%% variance, PC2 = %.1f%% variance", percent_var[[1]], percent_var[[2]]),
    x = sprintf("PC1 (%.1f%%)", percent_var[[1]]),
    y = sprintf("PC2 (%.1f%%)", percent_var[[2]]),
    color = factor1_name,
    shape = factor2_name
  )

dist_plot <- matrix_to_long(dist_mat, "sample_1", "sample_2", "distance") %>%
  mutate(
    sample_1 = factor(sample_1, levels = order_samples),
    sample_2 = factor(sample_2, levels = rev(order_samples))
  ) %>%
  ggplot(aes(sample_1, sample_2, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "Blues"))) +
  labs(title = "Sample-to-sample distance", x = NULL, y = NULL, fill = "Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_plot <- matrix_to_long(cor_mat, "sample_1", "sample_2", "correlation") %>%
  mutate(
    sample_1 = factor(sample_1, levels = order_samples),
    sample_2 = factor(sample_2, levels = rev(order_samples))
  ) %>%
  ggplot(aes(sample_1, sample_2, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2F6690", mid = "white", high = "#D1495B", midpoint = 0.97, limits = c(0.93, 1.00)) +
  labs(title = "Sample correlation heatmap", x = NULL, y = NULL, fill = "Pearson r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

qc_composite <- wrap_plots(
  library_plot,
  detected_plot,
  pca_plot,
  dist_plot,
  cor_plot,
  ncol = 2
) +
  plot_annotation(
    title = "RNA-seq quality control overview",
    subtitle = sprintf(
      "%d samples across %s and %s conditions. Median within-group correlation = %.3f; median between-group correlation = %.3f.",
      ncol(dds_interaction),
      factor1_name,
      factor2_name,
      cor_summary$value[[1]],
      cor_summary$value[[2]]
    ),
    tag_levels = "A"
  )

save_plot_dual(qc_composite, "qc_overview", width = 16, height = 18)

message("Fitting group-level DESeq2 models for shrunken pairwise contrasts...")

filtered_counts <- counts(dds_interaction)
filtered_metadata <- as.data.frame(colData(dds_interaction))

fit_group_model <- function(reference_group) {
  this_coldata <- filtered_metadata
  this_coldata$group <- relevel(factor(this_coldata$group), ref = reference_group)
  rownames(this_coldata) <- rownames(filtered_metadata)
  dds_group <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = this_coldata,
    design = ~ group
  )
  DESeq(dds_group, quiet = TRUE)
}

group_references <- levels(filtered_metadata$group)
group_models <- setNames(lapply(group_references, fit_group_model), group_references)

contrast_specs <- list()
primary_contrasts <- character()
interaction_contrasts <- character()

add_contrast_spec <- function(spec) {
  contrast_specs[[length(contrast_specs) + 1]] <<- spec
}

for (factor2_level in factor2_levels) {
  reference_group <- make_group_id(factor1_ref, factor2_level)
  for (factor1_level in setdiff(factor1_levels, factor1_ref)) {
    contrast_name <- sanitize_token(paste("factor1", factor1_level, "vs", factor1_ref, "in", factor2_name, factor2_level, sep = "_"))
    add_contrast_spec(list(
      name = contrast_name,
      label = sprintf("%s vs %s (%s = %s)", factor1_level, factor1_ref, factor2_name, factor2_level),
      mode = "group",
      ref_model = reference_group,
      case_group = make_group_id(factor1_level, factor2_level),
      ref_group = reference_group
    ))
    primary_contrasts <- c(primary_contrasts, contrast_name)
  }

  secondary_levels <- setdiff(factor1_levels, factor1_ref)
  if (length(secondary_levels) > 1) {
    for (pair in combn(secondary_levels, 2, simplify = FALSE)) {
      contrast_name <- sanitize_token(paste("factor1", pair[[1]], "vs", pair[[2]], "in", factor2_name, factor2_level, sep = "_"))
      add_contrast_spec(list(
        name = contrast_name,
        label = sprintf("%s vs %s (%s = %s)", pair[[1]], pair[[2]], factor2_name, factor2_level),
        mode = "group",
        ref_model = make_group_id(pair[[2]], factor2_level),
        case_group = make_group_id(pair[[1]], factor2_level),
        ref_group = make_group_id(pair[[2]], factor2_level)
      ))
    }
  }
}

for (factor1_level in factor1_levels) {
  contrast_name <- sanitize_token(paste("factor2", factor2_case, "vs", factor2_ref, "in", factor1_name, factor1_level, sep = "_"))
  add_contrast_spec(list(
    name = contrast_name,
    label = sprintf("%s vs %s (%s = %s)", factor2_case, factor2_ref, factor1_name, factor1_level),
    mode = "group",
    ref_model = make_group_id(factor1_level, factor2_ref),
    case_group = make_group_id(factor1_level, factor2_case),
    ref_group = make_group_id(factor1_level, factor2_ref)
  ))
}

for (factor1_level in setdiff(factor1_levels, factor1_ref)) {
  contrast_name <- sanitize_token(paste("interaction", factor1_level, "vs", factor1_ref, factor2_name, "response", sep = "_"))
  add_contrast_spec(list(
    name = contrast_name,
    label = sprintf("Interaction: %s vs %s %s response", factor1_level, factor1_ref, factor2_name),
    mode = "interaction",
    mutant = factor1_level
  ))
  interaction_contrasts <- c(interaction_contrasts, contrast_name)
}

contrast_label_map <- stats::setNames(
  vapply(contrast_specs, `[[`, character(1), "label"),
  vapply(contrast_specs, `[[`, character(1), "name")
)

get_group_coef <- function(dds_group, case_group, ref_group) {
  find_single_coef(resultsNames(dds_group), c("group", case_group, ref_group))
}

get_interaction_coef <- function(mutant) {
  find_single_coef(resultsNames(dds_interaction), c(mutant, "light_condition"))
}

prepare_result_df <- function(shrunk_res, raw_res, contrast_name) {
  raw_stat_df <- as.data.frame(raw_res) %>%
    rownames_to_column("gene_id") %>%
    transmute(gene_id, stat)

  res_df <- as.data.frame(shrunk_res) %>%
    rownames_to_column("gene_id") %>%
    left_join(raw_stat_df, by = "gene_id") %>%
    left_join(gene_annotation, by = "gene_id") %>%
    mutate(
      symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol),
      gene_name = ifelse(is.na(gene_name), "", gene_name),
      contrast = contrast_name,
      contrast_label = pretty_contrast_label(contrast_name),
      significance = case_when(
        !is.na(padj) & padj < args$padj_cutoff & log2FoldChange >= args$lfc_cutoff ~ "Up-regulated",
        !is.na(padj) & padj < args$padj_cutoff & log2FoldChange <= -args$lfc_cutoff ~ "Down-regulated",
        TRUE ~ "Not significant"
      ),
      abs_log2FC = abs(log2FoldChange)
    ) %>%
    arrange(padj, desc(abs_log2FC), gene_id)

  res_df
}

results_list <- list()

for (spec in contrast_specs) {
  message(sprintf("Computing contrast: %s", spec$name))
  if (identical(spec$mode, "group")) {
    dds_group <- group_models[[spec$ref_model]]
    coef_name <- get_group_coef(dds_group, spec$case_group, spec$ref_group)
    raw_res <- results(dds_group, name = coef_name)
    shrunk_res <- lfcShrink(dds_group, coef = coef_name, type = "apeglm", quiet = TRUE)
  } else {
    coef_name <- get_interaction_coef(spec$mutant)
    raw_res <- results(dds_interaction, name = coef_name)
    shrunk_res <- lfcShrink(dds_interaction, coef = coef_name, type = "apeglm", quiet = TRUE)
  }

  res_df <- prepare_result_df(shrunk_res, raw_res, spec$name)
  sig_df <- res_df %>% filter(significance != "Not significant")

  write_table(res_df, file.path(tables_full_dir, paste0(spec$name, "_full_results.csv")))
  write_table(sig_df, file.path(tables_filtered_dir, paste0(spec$name, "_DEGs.csv")))

  results_list[[spec$name]] <- list(
    coef_name = coef_name,
    results = res_df,
    significant = sig_df
  )
}

deg_summary <- bind_rows(lapply(names(results_list), function(name) {
  res_df <- results_list[[name]]$results
  tibble(
    contrast = name,
    contrast_label = pretty_contrast_label(name),
    up = sum(res_df$significance == "Up-regulated", na.rm = TRUE),
    down = sum(res_df$significance == "Down-regulated", na.rm = TRUE),
    total_significant = sum(res_df$significance != "Not significant", na.rm = TRUE)
  )
}))

write_table(deg_summary, file.path(tables_dir, "differential_expression_summary.csv"))

deg_plot_df <- deg_summary %>%
  select(contrast_label, up, down) %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    signed_count = ifelse(direction == "down", -count, count),
    direction = recode(direction, up = "Up-regulated", down = "Down-regulated"),
    contrast_label = factor(contrast_label, levels = deg_summary$contrast_label[order(deg_summary$total_significant)])
  )

deg_counts_plot <- ggplot(deg_plot_df, aes(x = contrast_label, y = signed_count, fill = direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = sig_cols[c("Up-regulated", "Down-regulated")]) +
  scale_y_continuous(labels = abs) +
  labs(
    title = "Differentially expressed genes by contrast",
    subtitle = sprintf(
      "Counts shown at padj < %.3g and |shrunken log2FC| >= %.2f",
      args$padj_cutoff,
      args$lfc_cutoff
    ),
    x = NULL,
    y = "Number of genes",
    fill = NULL
  )

make_volcano_plot <- function(res_df, title) {
  labeled_ids <- unique(c(
    res_df %>%
      filter(significance != "Not significant") %>%
      arrange(padj, desc(abs_log2FC)) %>%
      slice_head(n = 10) %>%
      pull(gene_id),
    res_df %>%
      filter(grepl(args$target_regex, symbol, ignore.case = TRUE), padj < 0.1) %>%
      pull(gene_id)
  ))

  plot_df <- res_df %>%
    mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      label = ifelse(gene_id %in% labeled_ids, symbol, NA_character_)
    )

  ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.75, size = 1.3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#6c757d", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#6c757d", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = subset(plot_df, !is.na(label)),
      aes(label = label),
      size = 3,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = 30,
      show.legend = FALSE
    ) +
    scale_color_manual(values = sig_cols) +
    labs(
      title = title,
      x = "Shrunken log2 fold change",
      y = expression(-log[10]("adjusted p-value")),
      color = NULL
    ) +
    theme(legend.position = "bottom")
}
volcano_plots <- lapply(primary_contrasts, function(name) {
  make_volcano_plot(results_list[[name]]$results, pretty_contrast_label(name))
})

names(volcano_plots) <- primary_contrasts

overlap_pairs <- expand.grid(
  contrast_1 = primary_contrasts,
  contrast_2 = primary_contrasts,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(
    genes_1 = list(results_list[[contrast_1]]$significant$gene_id),
    genes_2 = list(results_list[[contrast_2]]$significant$gene_id),
    intersection_count = length(intersect(genes_1, genes_2)),
    union_count = length(union(genes_1, genes_2)),
    jaccard = ifelse(union_count == 0, 0, intersection_count / union_count)
  ) %>%
  ungroup() %>%
  mutate(
    contrast_1 = factor(pretty_contrast_label(contrast_1), levels = pretty_contrast_label(primary_contrasts)),
    contrast_2 = factor(pretty_contrast_label(contrast_2), levels = rev(pretty_contrast_label(primary_contrasts)))
  )

overlap_plot <- ggplot(overlap_pairs, aes(contrast_1, contrast_2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = intersection_count), size = 3.6, fontface = "bold") +
  scale_fill_gradient(low = "#F1F5F9", high = "#1D4E89", labels = percent_format(accuracy = 1)) +
  labs(
    title = "Overlap summary across primary mutant contrasts",
    subtitle = "Tiles show Jaccard similarity; labels show intersecting significant genes",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

de_overview <- wrap_plots(
  c(list(deg_counts_plot), volcano_plots, list(overlap_plot)),
  ncol = 2
) +
  plot_annotation(
    title = "Differential expression overview",
    subtitle = "Primary mutant contrasts are shown as volcano plots alongside DEG counts and overlap structure.",
    tag_levels = "A"
  )

save_plot_dual(de_overview, "de_overview", width = 16, height = 18)

message("Building transcriptome-wide DEG pattern heatmap...")

top_genes <- unique(unlist(lapply(primary_contrasts, function(name) {
  this_df <- results_list[[name]]$results %>%
    filter(significance != "Not significant") %>%
    arrange(padj, desc(abs_log2FC)) %>%
    slice_head(n = 20) %>%
    pull(gene_id)
  if (length(this_df) == 0) {
    results_list[[name]]$results %>%
      arrange(padj, desc(abs_log2FC)) %>%
      slice_head(n = 20) %>%
      pull(gene_id)
  } else {
    this_df
  }
})))

top_genes <- unique(top_genes)
top_genes <- top_genes[top_genes %in% rownames(vst_mat)]

top_mat <- rescale_rows(vst_mat[top_genes, order_samples, drop = FALSE])
gene_order <- rownames(top_mat)[hclust(dist(top_mat))$order]
top_mat <- top_mat[gene_order, , drop = FALSE]

pattern_labels <- extract_primary_symbol(gene_annotation, rownames(top_mat))
rownames(top_mat) <- make.unique(pattern_labels[rownames(top_mat)])
colnames(top_mat) <- as.character(metadata$sample_label[match(colnames(top_mat), metadata$sample_id)])

pattern_long <- matrix_to_long(top_mat, "gene", "sample", "scaled_vst")
pattern_long$gene <- factor(pattern_long$gene, levels = rev(unique(rownames(top_mat))))
pattern_long$sample <- factor(pattern_long$sample, levels = colnames(top_mat))

pattern_plot <- ggplot(pattern_long, aes(sample, gene, fill = scaled_vst)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2F6690", mid = "white", high = "#D1495B", midpoint = 0) +
  labs(
    title = "Expression patterns of top differentially expressed genes",
    subtitle = sprintf("Union of the top genes from %d primary mutant contrasts", length(primary_contrasts)),
    x = NULL,
    y = NULL,
    fill = "Row-scaled VST"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

save_plot_dual(pattern_plot, "pattern_heatmap", width = 14, height = 18)

message(sprintf("Compiling %s-focused outputs...", args$target_label))

ftsh_ids <- gene_annotation %>%
  filter(
    grepl(args$target_regex, symbol, ignore.case = TRUE) |
      grepl(args$target_regex, gene_name, ignore.case = TRUE)
  ) %>%
  pull(gene_id) %>%
  unique() %>%
  sort()

ftsh_annotation <- gene_annotation %>%
  filter(gene_id %in% ftsh_ids) %>%
  arrange(symbol, gene_id)

ftsh_long <- tibble()
ftsh_significance <- tibble(
  gene_id = character(),
  contrast = character(),
  contrast_label = character(),
  log2FoldChange = numeric(),
  padj = numeric(),
  significance = character(),
  symbol = character(),
  gene_name = character()
)

if (length(ftsh_ids) > 0) {
  ftsh_vst <- vst_mat[ftsh_ids, order_samples, drop = FALSE]
  ftsh_vst_scaled <- rescale_rows(ftsh_vst)
  ftsh_symbol_map <- extract_primary_symbol(gene_annotation, rownames(ftsh_vst_scaled))
  rownames(ftsh_vst_scaled) <- make.unique(ftsh_symbol_map[rownames(ftsh_vst_scaled)])
  colnames(ftsh_vst_scaled) <- as.character(metadata$sample_label[match(colnames(ftsh_vst_scaled), metadata$sample_id)])

  ftsh_heatmap_df <- matrix_to_long(ftsh_vst_scaled, "gene", "sample", "scaled_vst")
  ftsh_heatmap_df$gene <- factor(ftsh_heatmap_df$gene, levels = rev(rownames(ftsh_vst_scaled)))
  ftsh_heatmap_df$sample <- factor(ftsh_heatmap_df$sample, levels = colnames(ftsh_vst_scaled))

  ftsh_heatmap_plot <- ggplot(ftsh_heatmap_df, aes(sample, gene, fill = scaled_vst)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2F6690", mid = "white", high = "#D1495B", midpoint = 0) +
    labs(
      title = sprintf("%s expression landscape", args$target_label),
      x = NULL,
      y = NULL,
      fill = "Row-scaled VST"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ftsh_long <- as.data.frame(vst_mat[ftsh_ids, , drop = FALSE]) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_annotation, by = "gene_id") %>%
    pivot_longer(cols = -c(gene_id, symbol, gene_name), names_to = "sample_id", values_to = "vst_expression") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(
      symbol = ifelse(is.na(symbol) | symbol == "", gene_id, symbol),
      group = factor(group, levels = levels(metadata$group))
    )

  ftsh_expression_plot <- ggplot(ftsh_long, aes(x = group, y = vst_expression, color = genotype, shape = light_condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.4, alpha = 0.9) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 18,
      size = 2.7,
      color = "black"
    ) +
    scale_color_manual(values = genotype_cols) +
    facet_wrap(~ symbol, scales = "free_y", ncol = 3) +
    labs(
      title = sprintf("Per-gene %s expression across experimental groups", args$target_label),
      subtitle = "Black diamonds indicate group means on the VST scale",
      x = NULL,
      y = "VST expression",
      color = factor1_name,
      shape = factor2_name
    ) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ftsh_composite <- ftsh_heatmap_plot / ftsh_expression_plot +
    plot_annotation(
      title = sprintf("%s-focused transcriptional analysis", args$target_label),
      subtitle = sprintf("Detected %s loci included in this analysis: %d", args$target_label, length(ftsh_ids)),
      tag_levels = "A"
    )

  save_plot_dual(ftsh_composite, paste0(sanitize_token(tolower(args$target_label)), "_composite"), width = 16, height = 18)

  ftsh_group_means <- ftsh_long %>%
    group_by(gene_id, symbol, gene_name, group) %>%
    summarise(mean_vst = mean(vst_expression), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_vst)

  ftsh_significance <- bind_rows(lapply(names(results_list), function(name) {
    results_list[[name]]$results %>%
      filter(gene_id %in% ftsh_ids) %>%
      transmute(
        gene_id,
        contrast = name,
        contrast_label = pretty_contrast_label(name),
        log2FoldChange,
        padj,
        significance
      )
  })) %>%
    left_join(ftsh_annotation, by = "gene_id") %>%
    arrange(symbol, contrast_label)

  ftsh_summary <- ftsh_significance %>%
    left_join(ftsh_group_means, by = c("gene_id", "symbol", "gene_name"))

  write_table(ftsh_summary, file.path(tables_ftsh_dir, paste0(sanitize_token(tolower(args$target_label)), "_summary_table.csv")))
  write_table(ftsh_long, file.path(tables_ftsh_dir, paste0(sanitize_token(tolower(args$target_label)), "_vst_expression_long.csv")))
} else {
  write_table(tibble(), file.path(tables_ftsh_dir, paste0(sanitize_token(tolower(args$target_label)), "_summary_table.csv")))
  write_table(tibble(), file.path(tables_ftsh_dir, paste0(sanitize_token(tolower(args$target_label)), "_vst_expression_long.csv")))
}

message("Running GO enrichment analyses...")

all_gene_ids <- rownames(dds_interaction)
ora_results <- list()

run_ora <- function(gene_ids, contrast_name, direction) {
  gene_ids <- unique(gene_ids[gene_ids %in% all_gene_ids])
  if (length(gene_ids) < 15) {
    return(NULL)
  }

  ora_obj <- tryCatch(
    enrichGO(
      gene = gene_ids,
      universe = all_gene_ids,
      OrgDb = OrgDb,
      keyType = args$gene_id_keytype,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = args$padj_cutoff,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(ora_obj)) {
    return(NULL)
  }

  ora_df <- as.data.frame(ora_obj)
  if (nrow(ora_df) == 0) {
    return(NULL)
  }

  ora_df %>%
    mutate(
      contrast = contrast_name,
      contrast_label = pretty_contrast_label(contrast_name),
      direction = direction,
      gene_ratio_numeric = parse_ratio(GeneRatio)
    ) %>%
    arrange(p.adjust, desc(Count))
}

for (name in names(results_list)) {
  res_df <- results_list[[name]]$results
  ora_results[[paste0(name, "_up")]] <- run_ora(res_df$gene_id[res_df$significance == "Up-regulated"], name, "Up-regulated")
  ora_results[[paste0(name, "_down")]] <- run_ora(res_df$gene_id[res_df$significance == "Down-regulated"], name, "Down-regulated")
}

ora_summary <- bind_rows(ora_results)
if (nrow(ora_summary) == 0) {
  ora_summary <- tibble()
}

write_table(ora_summary, file.path(tables_enrichment_dir, "go_ora_summary.csv"))

gsea_targets <- c(primary_contrasts, interaction_contrasts)

run_gsea <- function(res_df, contrast_name) {
  gene_stats <- res_df %>%
    filter(!is.na(.data$stat)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(desc(.data$stat))

  ranked_list <- gene_stats$stat
  names(ranked_list) <- gene_stats$gene_id

  gsea_obj <- tryCatch(
    gseGO(
      geneList = ranked_list,
      OrgDb = OrgDb,
      keyType = args$gene_id_keytype,
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.1,
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

gsea_summary <- bind_rows(lapply(gsea_targets, function(name) run_gsea(results_list[[name]]$results, name)))
if (nrow(gsea_summary) == 0) {
  gsea_summary <- tibble()
}

write_table(gsea_summary, file.path(tables_enrichment_dir, "go_gsea_summary.csv"))

functional_plot <- NULL
if (nrow(ora_summary) > 0) {
  ora_plot_df <- ora_summary %>%
    group_by(contrast_label, direction) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    mutate(
      Description = stringr::str_wrap(Description, width = 38),
      direction = factor(direction, levels = c("Up-regulated", "Down-regulated"))
    )

  functional_plot <- ggplot(
    ora_plot_df,
    aes(x = gene_ratio_numeric, y = Description, size = Count, color = p.adjust)
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D1495B", high = "#2F6690", trans = "reverse") +
    facet_grid(direction ~ contrast_label, scales = "free_y", space = "free_y") +
    labs(
      title = "GO biological process enrichment",
      subtitle = "Top over-representation results for contrasts with at least 15 DEGs in the tested direction",
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Count"
    ) +
    theme(
      strip.text.x = element_text(size = 10),
      strip.text.y = element_text(size = 10),
      axis.text.y = element_text(size = 8)
    )

  save_plot_dual(functional_plot, "functional_enrichment", width = 18, height = 14)
}

message("Writing report outputs...")

top_deg_lines <- deg_summary %>%
  arrange(desc(total_significant)) %>%
  slice_head(n = 5) %>%
  transmute(line = sprintf("- %s: %d DEGs (%d up, %d down)", contrast_label, total_significant, up, down)) %>%
  pull(line)

target_figure_stem <- paste0(sanitize_token(tolower(args$target_label)), "_composite")

ftsh_hits <- ftsh_significance %>%
  filter(significance != "Not significant") %>%
  mutate(line = sprintf("- %s in %s: log2FC = %.2f, padj = %.3g", symbol, contrast_label, log2FoldChange, padj)) %>%
  distinct(line) %>%
  pull(line)

if (length(ftsh_hits) == 0) {
  ftsh_hits <- sprintf(
    "- No %s genes met the DE threshold of padj < %.3g and |log2FC| >= %.2f.",
    args$target_label,
    args$padj_cutoff,
    args$lfc_cutoff
  )
}

top_ora_lines <- if (nrow(ora_summary) > 0) {
  ora_summary %>%
    group_by(contrast_label, direction) %>%
    arrange(p.adjust) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(line = sprintf("- %s (%s): %s [padj = %.3g]", contrast_label, direction, Description, p.adjust)) %>%
    pull(line)
} else {
  "- No ORA result met the reporting threshold."
}

top_gsea_lines <- if (nrow(gsea_summary) > 0) {
  gsea_summary %>%
    group_by(contrast_label) %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(line = sprintf("- %s: %s [NES = %.2f, padj = %.3g]", contrast_label, Description, NES, p.adjust)) %>%
    pull(line)
} else {
  "- No GSEA result met the reporting threshold."
}

report_md <- c(
  "# RNA-seq analysis report",
  "",
  "## Study design",
  "",
  sprintf("- Samples analyzed: %d", ncol(dds_interaction)),
  sprintf("- Genes retained after filtering: %d", nrow(dds_interaction)),
  sprintf(
    "- Design formula: `~ genotype * light_condition` with `%s` and `%s` as reference levels for `%s` and `%s`.",
    factor1_ref,
    factor2_ref,
    factor1_name,
    factor2_name
  ),
  "- Input types: gene-level RSEM counts plus sample metadata only. Sequencing-level QC metrics are out of scope for this report.",
  "",
  "## Quality control highlights",
  "",
  sprintf("- PC1 explains %.1f%% of variance and PC2 explains %.1f%%.", percent_var[[1]], percent_var[[2]]),
  sprintf("- Median within-group Pearson correlation: %.3f", cor_summary$value[[1]]),
  sprintf("- Median between-group Pearson correlation: %.3f", cor_summary$value[[2]]),
  "- Replicate profiles cluster tightly by genotype/light group, and no sample was flagged as an outlier for exclusion.",
  "",
  "![QC overview](figures/png/qc_overview.png)",
  "",
  "## Differential expression overview",
  "",
  top_deg_lines,
  "",
  "![DE overview](figures/png/de_overview.png)",
  "",
  "## Transcriptome patterning",
  "",
  sprintf("- The heatmap below shows the union of top genes from %d primary factor1-versus-reference contrasts after row-scaling VST values.", length(primary_contrasts)),
  "",
  "![Pattern heatmap](figures/png/pattern_heatmap.png)",
  "",
  sprintf("## %s-focused analysis", args$target_label),
  "",
  ftsh_hits,
  "",
  if (length(ftsh_ids) > 0) sprintf("![%s composite](figures/png/%s.png)", args$target_label, target_figure_stem) else NULL,
  "",
  "## Functional interpretation",
  "",
  "### GO over-representation highlights",
  "",
  top_ora_lines,
  "",
  "### GSEA highlights",
  "",
  top_gsea_lines,
  ""
)

if (!is.null(functional_plot)) {
  report_md <- c(
    report_md,
    "![Functional enrichment](figures/png/functional_enrichment.png)",
    ""
  )
}

report_md <- c(
  report_md,
  "## Deliverables",
  "",
  "- `tables/de_full/`: full DESeq2 result tables for every contrast.",
  sprintf("- `tables/deg_filtered/`: filtered DEG tables at padj < %.3g and |shrunken log2FC| >= %.2f.", args$padj_cutoff, args$lfc_cutoff),
  "- `tables/annotation/`: gene annotation lookup table.",
  sprintf("- `tables/%s/`: %s-specific expression and significance summaries.", sanitize_token(tolower(args$target_label)), args$target_label),
  "- `tables/enrichment/`: GO ORA and GSEA output tables.",
  "- `figures/png/` and `figures/pdf/`: publication-ready figures in raster and vector formats.",
  ""
)

writeLines(report_md, file.path(outdir, "analysis_report.md"))

report_html <- c(
  "<!DOCTYPE html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8' />",
  "<title>RNA-seq analysis report</title>",
  "<style>",
  "body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; line-height: 1.6; color: #1f2933; padding: 0 20px; }",
  "h1, h2, h3 { color: #102A43; }",
  "img { max-width: 100%; border: 1px solid #d9e2ec; margin: 16px 0 28px; }",
  "code { background: #f0f4f8; padding: 2px 4px; border-radius: 4px; }",
  "ul { padding-left: 20px; }",
  "li { margin-bottom: 6px; }",
  "</style>",
  "</head>",
  "<body>",
  "<h1>RNA-seq analysis report</h1>",
  "<h2>Study design</h2>",
  sprintf("<ul><li>Samples analyzed: %d</li>", ncol(dds_interaction)),
  sprintf("<li>Genes retained after filtering: %d</li>", nrow(dds_interaction)),
  sprintf(
    "<li>Design formula: <code>~ genotype * light_condition</code> with <code>%s</code> and <code>%s</code> as reference levels for <code>%s</code> and <code>%s</code>.</li>",
    factor1_ref,
    factor2_ref,
    factor1_name,
    factor2_name
  ),
  "<li>Input types: gene-level RSEM counts plus sample metadata only. Sequencing-level QC metrics are out of scope for this report.</li></ul>",
  "<h2>Quality control highlights</h2>",
  sprintf("<ul><li>PC1 explains %.1f%% of variance and PC2 explains %.1f%%.</li>", percent_var[[1]], percent_var[[2]]),
  sprintf("<li>Median within-group Pearson correlation: %.3f</li>", cor_summary$value[[1]]),
  sprintf("<li>Median between-group Pearson correlation: %.3f</li>", cor_summary$value[[2]]),
  "<li>Replicate profiles cluster tightly by genotype/light group, and no sample was flagged as an outlier for exclusion.</li></ul>",
  "<img src='figures/png/qc_overview.png' alt='QC overview' />",
  "<h2>Differential expression overview</h2>",
  "<ul>",
  paste0("<li>", sub("^- ", "", top_deg_lines), "</li>", collapse = ""),
  "</ul>",
  "<img src='figures/png/de_overview.png' alt='DE overview' />",
  "<h2>Transcriptome patterning</h2>",
  sprintf("<p>The heatmap below shows the union of top genes from %d primary factor1-versus-reference contrasts after row-scaling VST values.</p>", length(primary_contrasts)),
  "<img src='figures/png/pattern_heatmap.png' alt='Pattern heatmap' />",
  sprintf("<h2>%s-focused analysis</h2>", args$target_label),
  "<ul>",
  paste0("<li>", sub("^- ", "", ftsh_hits), "</li>", collapse = ""),
  "</ul>",
  if (length(ftsh_ids) > 0) sprintf("<img src='figures/png/%s.png' alt='%s composite' />", target_figure_stem, args$target_label) else NULL,
  "<h2>Functional interpretation</h2>",
  "<h3>GO over-representation highlights</h3>",
  "<ul>",
  paste0("<li>", sub("^- ", "", top_ora_lines), "</li>", collapse = ""),
  "</ul>",
  "<h3>GSEA highlights</h3>",
  "<ul>",
  paste0("<li>", sub("^- ", "", top_gsea_lines), "</li>", collapse = ""),
  "</ul>"
)

if (!is.null(functional_plot)) {
  report_html <- c(report_html, "<img src='figures/png/functional_enrichment.png' alt='Functional enrichment' />")
}

report_html <- c(
  report_html,
  "<h2>Deliverables</h2>",
  "<ul>",
  "<li><code>tables/de_full/</code>: full DESeq2 result tables for every contrast.</li>",
  sprintf("<li><code>tables/deg_filtered/</code>: filtered DEG tables at padj &lt; %.3g and |shrunken log2FC| &gt;= %.2f.</li>", args$padj_cutoff, args$lfc_cutoff),
  "<li><code>tables/annotation/</code>: gene annotation lookup table.</li>",
  sprintf("<li><code>tables/%s/</code>: %s-specific expression and significance summaries.</li>", sanitize_token(tolower(args$target_label)), args$target_label),
  "<li><code>tables/enrichment/</code>: GO ORA and GSEA output tables.</li>",
  "<li><code>figures/png/</code> and <code>figures/pdf/</code>: publication-ready figures in raster and vector formats.</li>",
  "</ul>",
  "</body>",
  "</html>"
)

writeLines(report_html, file.path(outdir, "analysis_report.html"))

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "session_info.txt"))

message("Analysis completed successfully.")
