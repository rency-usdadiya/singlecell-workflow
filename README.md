# singlecell-workflow

This repository contains a only workflow of single-cell RNA-seq preprocessing pipeline using Seurat, scDblFinder and SingleR.

Important: This repo intentionally does not include any raw input data or processed result files. Keep data local or host it externally and point the pipeline to your files.

What is included
- `script.R` — the complete analysis workflow (read, QC, doublet detection, clustering, annotation).

# Design principle
- Workflow-only: code and configuration only; no raw data or large outputs stored here.
- Reproducible: the script is deterministic and can be run on any machine with dependencies installed.
- Safe defaults: by default `save_results <- FALSE` to avoid accidentally committing large files.

# How to run
1. Place your input files (matrix, barcodes, genes) somewhere local or on a data server you control.
2. Provide file paths via environment variables.
