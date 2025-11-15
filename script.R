# Single-cell preprocessing pipeline using Seurat, scDblFinder, SingleR
# Author: Rency Usdadiya
# Date: 2025-11-15

# DESCRIPTION ---------------------------------------------------------------
# This script implements the analysis workflow but does NOT bundle data or
# write large result files by default. Set `save_results <- TRUE` below to
# enable writing outputs locally (these will not be pushed to the repo).

# 0. Options ---------------------------------------------------------------
save_results <- FALSE  # keep FALSE for workflow-only repository
outdir <- "results"   # used only if save_results == TRUE

# 1. Libraries ---------------------------------------------------------------
install.packages(c("Seurat","dplyr","ggplot2","Matrix"))
install.packages("BiocManager")
BiocManager::install(c("scDblFinder","SingleCellExperiment","BiocParallel","celldex","SingleR"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(celldex)
library(SingleR)

# 2. User-provided inputs --------------------------------------------------
# Provide your input files (not included in this repo). Example:
# data_dir <- "path/to/your/data"
# matrix_file   <- file.path(data_dir, "matrix.mtx.gz")
# barcodes_file <- file.path(data_dir, "barcodes.tsv.gz")
# genes_file    <- file.path(data_dir, "genes.tsv.gz")

# Example of reading from environment variables (recommended):
matrix_file   <- Sys.getenv("SC_MATRIX_FILE")   # full path to matrix.mtx(.gz)
barcodes_file <- Sys.getenv("SC_BARCODES_FILE")
genes_file    <- Sys.getenv("SC_GENES_FILE")

if (matrix_file == "" || barcodes_file == "" || genes_file == "") {
  stop("Input files not provided. Set SC_MATRIX_FILE, SC_BARCODES_FILE and SC_GENES_FILE environment variables, or edit the script's input paths.")
}

# QC thresholds (change as per the requirement) -----------------------------------------
min_features <- 500
max_features <- 6000
max_mt_pct   <- 40
min_cells_feature <- 3
min_features_cell  <- 200

# 3. Read input -------------------------------------------------------------
cat("Reading matrix and annotation files...
")
sc_matrix <- Matrix::readMM(gzfile(matrix_file))
barcodes  <- read.table(gzfile(barcodes_file), header = FALSE, stringsAsFactors = FALSE)
genes     <- read.table(gzfile(genes_file), header = FALSE, stringsAsFactors = FALSE)

# 4. Clean genes & set row/col names ----------------------------------------
genes <- genes[, 1:2]
colnames(genes) <- c("ensembl_id", "symbol")
genes$ensembl_id <- trimws(as.character(genes$ensembl_id))
genes$symbol      <- trimws(as.character(genes$symbol))
# Remove version numbers from Ensembl IDs
genes$ensembl_id <- sub("\..*$", "", genes$ensembl_id)

# Sanity check
if (nrow(sc_matrix) != nrow(genes)) {
  stop(sprintf("Row mismatch: matrix rows = %d, genes rows = %d. Check inputs.", nrow(sc_matrix), nrow(genes)))
}

unique_symbols <- make.unique(ifelse(genes$symbol == "" | is.na(genes$symbol), genes$ensembl_id, genes$symbol))
rownames(sc_matrix) <- unique_symbols
colnames(sc_matrix) <- barcodes$V1

# 5. Create Seurat object ---------------------------------------------------
seurat_obj <- CreateSeuratObject(
  counts = sc_matrix,
  project = "Ovarian",
  min.cells = min_cells_feature,
  min.features = min_features_cell
)

# 6. QC ---------------------------------------------------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= min_features &
           nFeature_RNA <= max_features &
           percent.mt <= max_mt_pct
)

# 7. Doublet detection (scDblFinder) ---------------------------------------
cat("Converting to SingleCellExperiment and running scDblFinder...
")
sce <- as.SingleCellExperiment(seurat_obj)
# Use a single worker in cross-platform workflows for safety
sce <- scDblFinder(sce, BPPARAM = SnowParam(workers = 1, type = "SOCK"))
seurat_obj$scDblFinder.class <- colData(sce)$scDblFinder.class
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# 8. Normalization, HVG, scaling, PCA --------------------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# 9. Clustering & embeddings ------------------------------------------------
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# 10. Marker discovery -------------------------------------------------------
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)

# 11. Cell-type annotation with SingleR ------------------------------------
cat("Loading HPCA ref and running SingleR...
")
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
sce_for_singleR <- as.SingleCellExperiment(seurat_obj)
singleR_out <- SingleR(test = sce_for_singleR, ref = hpca.ref,
                       labels = hpca.ref$label.main,
                       clusters = seurat_obj$seurat_clusters)
cluster2label <- singleR_out$labels
names(cluster2label) <- rownames(singleR_out)
cell_cluster <- as.character(seurat_obj$seurat_clusters)
seurat_obj$SingleR_labels <- cluster2label[cell_cluster]

# 12. Plots (display only) -------------------------------------------------
print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP clusters"))
print(DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_labels", label = TRUE, repel = TRUE) + ggtitle("Cell Type Annotation (SingleR)"))

# 13. Optional saving -------
if (isTRUE(save_results)) {
  if (!dir.exists(outdir)) dir.create(outdir)
  saveRDS(seurat_obj, file = file.path(outdir, "seurat_obj_filtered.rds"))
  write.csv(cluster_markers, file = file.path(outdir, "cluster_markers.csv"), row.names = FALSE)
  cat("Results written to:", normalizePath(outdir), "
")
} else {
  cat("save_results is FALSE — no large outputs were written. Set save_results <- TRUE to save locally.
")
}
