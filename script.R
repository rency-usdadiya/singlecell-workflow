# Single-cell preprocessing pipeline using Seurat, scDblFinder, SingleR
# Author: Rency Usdadiya
# Date: 2025-11-15

# 0. Libraries ---------------------------------------------------------------
# install.packages(c("Seurat","dplyr","ggplot2","Matrix"))
# install.packages("BiocManager")
# BiocManager::install(c("scDblFinder","SingleCellExperiment","BiocParallel","celldex","SingleR"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(celldex)
library(SingleR)

# 1. Define Path -----------------------------------------------------
data_dir <- "C:/Users/RENCY/Desktop/data1/GSE184880_1"
matrix_file   <- file.path(data_dir, "GSM5599225_cancer1.matrix.mtx.gz")
barcodes_file <- file.path(data_dir, "GSM5599225_cancer1.barcodes.tsv.gz")
genes_file    <- file.path(data_dir, "GSM5599225_cancer1.genes.tsv.gz")

# QC thresholds (change as per the requirement)
min_features <- 500
max_features <- 6000
max_mt_pct   <- 40
min_cells_feature <- 3
min_features_cell  <- 200

# 2. Read input -------------------------------------------------------------
cat("Reading matrix and annotation files...\n")
sc_matrix <- Matrix::readMM(gzfile(matrix_file))
barcodes  <- read.table(gzfile(barcodes_file), header = FALSE, stringsAsFactors = FALSE)
genes     <- read.table(gzfile(genes_file), header = FALSE, stringsAsFactors = FALSE)

# 3. Clean genes & set row/col names ----------------------------------------
genes <- genes[, 1:2]
colnames(genes) <- c("ensembl_id", "symbol")
genes$ensembl_id <- trimws(as.character(genes$ensembl_id))
genes$symbol      <- trimws(as.character(genes$symbol))

# Remove version numbers from Ensembl IDs (e.g., ENSG000001.5 -> ENSG000001)
genes$ensembl_id <- sub("\\..*$", "", genes$ensembl_id)

# Sanity check: rows must match
if (nrow(sc_matrix) != nrow(genes)) {
  stop(sprintf("Row mismatch: matrix rows = %d, genes rows = %d. Check inputs.", nrow(sc_matrix), nrow(genes)))
}

# Make gene symbols unique and assign
unique_symbols <- make.unique(ifelse(genes$symbol == "" | is.na(genes$symbol), genes$ensembl_id, genes$symbol))
rownames(sc_matrix) <- unique_symbols
colnames(sc_matrix) <- barcodes$V1

# 4. Create Seurat object ---------------------------------------------------
seurat_obj <- CreateSeuratObject(
  counts = sc_matrix,
  project = "Ovarian",
  min.cells = min_cells_feature,
  min.features = min_features_cell
)

# 5. QC ---------------------------------------------------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt <= max_mt_pct
)

# 6. Doublet detection (scDblFinder) ---------------------------------------
cat("Converting to SingleCellExperiment and running scDblFinder...\n")
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce, BPPARAM = SnowParam(workers = 1, type = "SOCK"))
seurat_obj$scDblFinder.class <- colData(sce)$scDblFinder.class

# Keep singlets only
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# 7. Normalization, HVG, scaling, PCA --------------------------------------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# 8. Clustering & embeddings ------------------------------------------------
ElbowPlot(seurat_obj)   # interactive inspection recommended
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# 9. Marker discovery -------------------------------------------------------
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.25)

# 10. Cell-type annotation with SingleR ------------------------------------
cat("Loading HPCA ref and running SingleR...\n")
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
sce_for_singleR <- as.SingleCellExperiment(seurat_obj)

# Run cluster-level SingleR annotation (faster & more robust)
singleR_out <- SingleR(test = sce_for_singleR, ref = hpca.ref,
                       labels = hpca.ref$label.main,
                       clusters = seurat_obj$seurat_clusters)

# Map SingleR cluster labels back to cells
cluster2label <- singleR_out$labels
names(cluster2label) <- rownames(singleR_out)  # cluster names -> labels
# Create per-cell label vector
cell_cluster <- as.character(seurat_obj$seurat_clusters)
seurat_obj$SingleR_labels <- cluster2label[cell_cluster]

# 11. Plots -----------------------------------------------------------------
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("UMAP clusters")
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_labels", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Annotation (SingleR)")

print(p1)
print(p2)

# 12. Save results ----------------------------------------------------------
outdir <- "results"
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(seurat_obj, file = file.path(outdir, "seurat_obj_filtered.rds"))
write.csv(cluster_markers, file = file.path(outdir, "cluster_markers.csv"), row.names = FALSE)

cat("Pipeline finished. Results in:", normalizePath(outdir), "\n")