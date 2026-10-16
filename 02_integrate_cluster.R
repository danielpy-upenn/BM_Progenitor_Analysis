###############################################################################
# 02_integrate_cluster.R
#
# What this script does
#   1. Normalizes (LogNormalize) and finds 2000 variable genes (vst) per run
#      layer -- the Seurat v5 equivalent of the per-sample NormalizeData /
#      FindVariableFeatures calls in preprocessing.R.
#   2. Scales, runs PCA, and makes an unintegrated UMAP (to show batch effects).
#   3. Integrates all runs with CCA:
#        IntegrateLayers(method = CCAIntegration)  (v5 version of v3 IntegrateData)
#   4. Clusters at resolutions 0.2-2.0 and runs UMAP with the paper's settings
#      (n.neighbors = 40, min.dist = 0.3).
#   5. Joins the RNA layers back together (needed by BoneMarrowMap and DE).
#   6. Saves UMAPs by cluster (one file per resolution), run, donor,
#      condition, sex and library type, plus cluster composition tables.
#
# Input  : Output/checkpoints/01_qc_filtered_merged.rds
# Outputs: Output/02_integration/, Output/checkpoints/02_integrated_clustered.rds
#
# Run    : Rscript 02_integrate_cluster.R
# Note   : CCA over 30 runs is the slowest, most memory-hungry step (hours).
###############################################################################

get_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
source(file.path(get_script_dir(), "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

UMAP_DIR    <- file.path(INT_DIR, "umap_by_metadata")
CLUSTER_DIR <- file.path(INT_DIR, "umap_by_resolution")
ensure_dirs(INT_DIR, UMAP_DIR, CLUSTER_DIR, LOG_DIR)
log_msg("=== 02_integrate_cluster.R started ===")
gc(reset = TRUE)
set.seed(SEED)

if (!file.exists(CKPT_01)) stop("Checkpoint from script 01 not found: ", CKPT_01)
obj <- readRDS(CKPT_01)
log_msg("Loaded %d genes x %d cells, %d runs", nrow(obj), ncol(obj), length(unique(obj$Run)))

###############################################################################
# Step 1. Normalize, variable features, scale, PCA
###############################################################################
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_VARIABLE_FEATURES)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = N_PCS_COMPUTE, seed.use = SEED, verbose = FALSE)
log_mem("after PCA")

p_elbow <- ElbowPlot(obj, ndims = N_PCS_COMPUTE, reduction = "pca") +
  geom_vline(xintercept = max(DIMS_USE), linetype = "dashed", colour = "#C0392B") +
  ggtitle(sprintf("Elbow plot (dashed = %d PCs used)", max(DIMS_USE))) + theme_classic()
save_plot(p_elbow, INT_DIR, "elbow_plot", width = 6, height = 4.5)

# Unintegrated UMAP, only to judge how much batch effect CCA removes
obj <- RunUMAP(obj, reduction = "pca", dims = DIMS_USE, n.neighbors = UMAP_N_NEIGHBORS,
               min.dist = UMAP_MIN_DIST, reduction.name = "umap.unintegrated",
               seed.use = SEED, verbose = FALSE)
p_unint <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "Run", raster = TRUE) +
  ggtitle("Before integration, colored by run") + theme(legend.text = element_text(size = 6))
save_plot(p_unint, INT_DIR, "umap_UNINTEGRATED_by_run", width = 10, height = 7)
p_unint_lt <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "LibraryType", raster = TRUE) +
  ggtitle("Before integration, colored by library type")
save_plot(p_unint_lt, INT_DIR, "umap_UNINTEGRATED_by_librarytype", width = 8, height = 6.5)

###############################################################################
# Step 2. CCA integration across runs
###############################################################################
cells_per_run <- table(obj$Run)
min_cells <- min(cells_per_run)
k_weight <- CCA_K_WEIGHT_DEFAULT
if (min_cells <= k_weight) {
  k_weight <- max(10, min_cells - 1)
  log_msg("Smallest run (%s) has %d cells; lowering k.weight from %d to %d",
          names(which.min(cells_per_run)), min_cells, CCA_K_WEIGHT_DEFAULT, k_weight)
}

reference_idx <- NULL
if (!is.null(INTEGRATION_REFERENCE_RUNS)) {
  layer_runs <- sub("^counts\\.", "", Layers(obj, search = "counts"))
  reference_idx <- match(INTEGRATION_REFERENCE_RUNS, layer_runs)
  if (anyNA(reference_idx)) stop("INTEGRATION_REFERENCE_RUNS contains runs not in the object")
  log_msg("Using reference-based integration with %d reference runs", length(reference_idx))
}

log_msg("Starting CCA integration across %d runs (k.weight = %d). This can take hours.",
        length(cells_per_run), k_weight)
t0 <- Sys.time()
obj <- IntegrateLayers(object = obj, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       dims = DIMS_USE, k.weight = k_weight,
                       reference = reference_idx, verbose = TRUE)
log_msg("Integration finished in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_mem("after integration")

###############################################################################
# Step 3. Neighbors, clusters (several resolutions), UMAP
###############################################################################
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = DIMS_USE, verbose = FALSE)
obj <- FindClusters(obj, resolution = CLUSTER_RESOLUTIONS, random.seed = SEED, verbose = FALSE)
res_cols <- paste0("RNA_snn_res.", CLUSTER_RESOLUTIONS)
missing_res <- setdiff(res_cols, colnames(obj@meta.data))
if (length(missing_res) > 0) {
  # Seurat drops trailing zeros in column names (e.g. "1" instead of "1.0")
  res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
}
log_msg("Cluster columns: %s", paste(res_cols, collapse = ", "))

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = DIMS_USE,
               n.neighbors = UMAP_N_NEIGHBORS, min.dist = UMAP_MIN_DIST,
               reduction.name = "umap", seed.use = SEED, verbose = FALSE)

# Rejoin layers: one counts + one data layer (split layers are only needed for
# integration). The integrated.cca reduction and the UMAP are kept.
obj <- JoinLayers(obj)
log_mem("after UMAP and JoinLayers")

###############################################################################
# Step 4. Plots
###############################################################################
# One labeled UMAP per resolution
for (rc in res_cols) {
  res <- res_label(rc)
  n_cl <- length(unique(obj[[rc, drop = TRUE]]))
  p <- DimPlot(obj, reduction = "umap", group.by = rc, label = TRUE, label.size = 3,
               repel = TRUE, raster = TRUE) + NoLegend() +
    ggtitle(sprintf("Resolution %s (%d clusters)", res, n_cl))
  save_plot(p, CLUSTER_DIR, sprintf("umap_clusters_res_%s", res), width = 8, height = 7.5)
}

meta_cols <- c("Run", "Donor", "Condition", "Sex", "LibraryType")
for (mc in meta_cols) {
  p <- DimPlot(obj, reduction = "umap", group.by = mc, raster = TRUE, shuffle = TRUE, seed = SEED) +
    ggtitle(sprintf("Integrated UMAP colored by %s", mc))
  if (mc == "Condition") p <- p + scale_colour_manual(values = CONDITION_COLORS)
  if (mc == "Run") p <- p + theme(legend.text = element_text(size = 6))
  save_plot(p, UMAP_DIR, paste0("umap_by_", mc), width = 9, height = 7)
}

# Split views: each donor / library type in its own panel (batch mixing check)
obj$Donor <- factor(obj$Donor, levels = donor_levels(obj$Donor))
p_split_cond <- DimPlot(obj, reduction = "umap", split.by = "Condition", group.by = "Condition",
                        raster = TRUE) + scale_colour_manual(values = CONDITION_COLORS) + NoLegend()
save_plot(p_split_cond, UMAP_DIR, "umap_split_by_Condition", width = 12, height = 6)
p_split_lt <- DimPlot(obj, reduction = "umap", split.by = "LibraryType", group.by = "LibraryType",
                      raster = TRUE, ncol = 2) + NoLegend()
save_plot(p_split_lt, UMAP_DIR, "umap_split_by_LibraryType", width = 11, height = 10)
p_split_donor <- DimPlot(obj, reduction = "umap", split.by = "Donor", group.by = "Condition",
                         raster = TRUE, ncol = 5) + scale_colour_manual(values = CONDITION_COLORS)
save_plot(p_split_donor, UMAP_DIR, "umap_split_by_Donor", width = 18, height = 11)
obj$Donor <- as.character(obj$Donor)

# QC metrics on the UMAP (helps spot low-quality clusters)
p_qc <- FeaturePlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    reduction = "umap", raster = TRUE, ncol = 3, order = FALSE)
save_plot(p_qc, UMAP_DIR, "umap_QC_metrics", width = 16, height = 5)

###############################################################################
# Step 5. Cluster composition tables (long format, all resolutions)
###############################################################################
comp <- lapply(res_cols, function(rc) {
  obj@meta.data %>%
    count(Cluster = .data[[rc]], Donor, Condition, LibraryType, name = "n_cells") %>%
    mutate(Resolution = res_label(rc), .before = 1)
}) %>% bind_rows()
write.csv(comp, file.path(INT_DIR, "cluster_composition_all_resolutions.csv"), row.names = FALSE)

cl_summary <- comp %>% group_by(Resolution) %>%
  summarise(n_clusters = n_distinct(Cluster), .groups = "drop")
write.csv(cl_summary, file.path(INT_DIR, "n_clusters_per_resolution.csv"), row.names = FALSE)
print(as.data.frame(cl_summary))

###############################################################################
# Save checkpoint, remove the 01 checkpoint
###############################################################################
saveRDS(obj, CKPT_02)
log_msg("Saved checkpoint %s", CKPT_02)
delete_checkpoint(CKPT_01)
log_mem("end of 02")
log_disk()
write_session_info("02_integrate_cluster")
log_msg("=== 02_integrate_cluster.R finished ===")
