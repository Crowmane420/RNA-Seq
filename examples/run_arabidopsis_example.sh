#!/usr/bin/env bash
set -euo pipefail

bash scripts/run_pipeline.sh \
  --counts=/path/to/rsem.merged.gene_counts.tsv \
  --metadata=/path/to/sample_metadata.csv \
  --outdir=/path/to/rnaseq_analysis \
  --libdir=/path/to/.Rlibs \
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
