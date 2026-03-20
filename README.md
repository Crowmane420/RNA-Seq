# RNA-seq DESeq2 Pipeline

A reusable two-factor RNA-seq analysis pipeline built around DESeq2, GO/KEGG enrichment, target-family tracking, WGCNA, regulator summaries, and exploratory GRN inference.

The pipeline was generalized from a validated Arabidopsis `factor1 x factor2` analysis, but it is designed to work with any organism and any count matrix plus metadata table that can be mapped into the expected schema.

## What It Produces

- QC summaries and publication-ready figures
- DESeq2 differential expression across:
  - all factor1-vs-reference contrasts within each factor2 level
  - all non-reference factor1 pairwise contrasts within each factor2 level
  - factor2 contrasts within each factor1 level
  - factor1-specific factor2 response interaction contrasts
- GO ORA and GSEA
- optional KEGG ORA and GSEA
- target-family expression summaries based on a regex
- regulator-focused summaries for transcription factors and kinases
- exploratory WGCNA and focused GRN outputs

## Current Scope

This repository currently supports:

- gene-level count matrices
- a two-factor experimental design
- exactly 2 levels for `factor2`
- any number of levels for `factor1`

It does not currently perform raw FASTQ processing, alignment, or transcript quantification. It assumes you already have a gene-level count matrix.

## Input Requirements

### Count matrix

- Tab-delimited text file
- One row per gene
- One column containing gene IDs
- Optional transcript column
- Remaining columns are sample counts

### Metadata

- CSV file
- One row per sample
- Must include:
  - sample ID column
  - factor1 column
  - factor2 column
- Optional:
  - replicate column
  - precomputed group column

If your metadata uses different column names, pass them via command-line arguments.

## Quick Start

Install dependencies:

```bash
Rscript scripts/install_dependencies.R --orgdb_package=org.At.tair.db
```

Run the full pipeline:

```bash
bash scripts/run_pipeline.sh \
  --counts=/path/to/counts.tsv \
  --metadata=/path/to/sample_metadata.csv \
  --outdir=/path/to/output_dir \
  --libdir=/path/to/Rlibs \
  --gene_id_col=gene_id \
  --transcript_id_col='transcript_id(s)' \
  --sample_id_col=sample_id \
  --factor1_col=genotype \
  --factor2_col=light_condition \
  --factor1_ref=WT \
  --factor2_ref=LL \
  --replicate_col=replicate \
  --group_col=group \
  --orgdb_package=org.At.tair.db \
  --gene_id_keytype=TAIR \
  --kegg_organism=ath \
  --target_regex='FTSH|FtsH' \
  --target_label=FTSH
```

## Main Arguments

- `--counts`: path to the gene-level count matrix
- `--metadata`: path to the sample metadata CSV
- `--outdir`: output directory
- `--libdir`: R package library directory
- `--gene_id_col`: gene ID column in the count matrix
- `--transcript_id_col`: transcript column in the count matrix; set to an empty string if absent
- `--sample_id_col`: sample ID column in the metadata
- `--factor1_col`: primary grouping factor
- `--factor2_col`: secondary grouping factor; must have exactly 2 levels
- `--factor1_ref`: reference level for factor1
- `--factor2_ref`: reference level for factor2
- `--replicate_col`: replicate column; set to an empty string if absent
- `--group_col`: optional precomputed group column
- `--orgdb_package`: organism annotation package, for example `org.At.tair.db` or `org.Hs.eg.db`
- `--gene_id_keytype`: key type understood by the selected OrgDb package
- `--kegg_organism`: KEGG organism code; set to an empty string to skip KEGG
- `--target_regex`: regex used to define the target gene family or pathway of interest
- `--target_label`: display label for the target gene family
- `--min_count`: count threshold for gene filtering
- `--min_samples`: sample threshold for gene filtering
- `--padj_cutoff`: adjusted p-value cutoff for DEG calls
- `--lfc_cutoff`: absolute shrunken log2 fold-change cutoff for DEG calls

## Output Structure

- `figures/png/` and `figures/pdf/`: publication-ready figures
- `tables/de_full/`: full DESeq2 results
- `tables/deg_filtered/`: filtered DEG tables
- `tables/annotation/`: annotation lookup tables
- `tables/enrichment/`: GO outputs
- `tables/kegg/`: KEGG outputs
- `tables/regulators/`: TF and kinase outputs
- `tables/wgcna/`: WGCNA outputs
- `tables/grn/`: focused GRN outputs
- `tables/<target_label>/`: target-family summaries
- `analysis_report.md` and `analysis_report.html`
- `advanced_regulatory_report.md` and `advanced_regulatory_report.html`

## Notes and Caveats

- The advanced extension is exploratory, especially WGCNA and GENIE3.
- KEGG support depends on the availability of a valid KEGG organism code and pathway mappings.
- Gene-family tracking depends on the quality of your annotation and the specificity of `--target_regex`.
- The pipeline is most appropriate for experiments with biologically meaningful replication and a clear two-factor design.
