###############################################################################
# 00_config.R
#
# Shared settings and helper functions for the FPD scRNA-seq pipeline.
# This file is NOT run on its own: every numbered script (01-04) sources it at
# the top, so paths and parameters live in exactly one place.
#
# To change a parameter (e.g. clustering resolutions, downsampling N),
# edit it here and rerun the affected script(s).
###############################################################################

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
# BASE_DIR can be overridden with an environment variable (used only for
# testing); on the EC2 node it defaults to /home/danielpy/data.
BASE_DIR   <- Sys.getenv("FPD_BASE_DIR", unset = "/home/danielpy/data")

SRR_DIR    <- file.path(BASE_DIR, "SRR_Folders")        # 30 cellranger folders
SRA_TABLE  <- file.path(BASE_DIR, "SraRunTable.csv")    # SRA run table
META_DIR   <- file.path(BASE_DIR, "metadata")           # small hand-edited tables
BMM_DIR    <- file.path(BASE_DIR, "BoneMarrowMap_ref")  # BoneMarrowMap reference files
OUT_DIR    <- file.path(BASE_DIR, "Output")

# Hand-edited input tables (all in META_DIR)
QC_THRESHOLDS_CSV     <- file.path(META_DIR, "qc_thresholds.csv")        # from preprocessing.R
RUN_LIBRARY_MAP_CSV   <- file.path(META_DIR, "run_library_map.csv")      # YOU fill this in
RUN_LIBRARY_TEMPLATE  <- file.path(META_DIR, "run_library_map_TEMPLATE.csv")
DOWNSAMPLE_OVERRIDES  <- file.path(META_DIR, "downsample_overrides.csv") # YOU fill this in
DOWNSAMPLE_TEMPLATE   <- file.path(META_DIR, "downsample_overrides_TEMPLATE.csv")

# BoneMarrowMap reference files (upload these yourself with the AWS CLI)
BMM_REF_RDS   <- file.path(BMM_DIR, "BoneMarrowMap_SymphonyReference.rds")
BMM_UWOT      <- file.path(BMM_DIR, "BoneMarrowMap_uwot_model.uwot")

# Output sub-folders
QC_DIR    <- file.path(OUT_DIR, "01_qc")
INT_DIR   <- file.path(OUT_DIR, "02_integration")
BMM_OUT   <- file.path(OUT_DIR, "03_bonemarrowmap")
DE_DIR    <- file.path(OUT_DIR, "04_de")
LOG_DIR   <- file.path(OUT_DIR, "logs")
CKPT_DIR  <- file.path(OUT_DIR, "checkpoints")

# Checkpoints (intermediate ones are deleted once the next stage succeeds)
CKPT_01   <- file.path(CKPT_DIR, "01_qc_filtered_merged.rds")
CKPT_02   <- file.path(CKPT_DIR, "02_integrated_clustered.rds")
FINAL_OBJ <- file.path(CKPT_DIR, "03_FPD_BM_annotated_FINAL.rds")   # kept

DELETE_INTERMEDIATE_CHECKPOINTS <- TRUE

# -----------------------------------------------------------------------------
# 2. Study design checks (script 01 stops if the run table disagrees)
# -----------------------------------------------------------------------------
EXPECTED_RUNS      <- sprintf("SRR%d", 29249442:29249471)   # 30 runs
EXPECTED_N_FPD     <- 10
EXPECTED_N_HEALTHY <- 4

# Known typos in submitted_subject_id: c(wrong = "right").
# FD19006 is a control (subject_is_affected == "No") and should read HD19006.
SUBJECT_ID_FIXES <- c(FD19006 = "HD19006")

# Allowed library types (MNC = unenriched; Enriched = "Spiked" = MNC + CD34+ + MSC)
LIBRARY_TYPES <- c("MNC", "Enriched", "CD11b+CD14+", "CD11b+CD14-")

# -----------------------------------------------------------------------------
# 3. QC parameters
# -----------------------------------------------------------------------------
MIN_CELLS_PER_GENE     <- 3      # CreateSeuratObject(min.cells)
MIN_FEATURES_PER_CELL  <- 200    # CreateSeuratObject(min.features)
MT_PATTERN             <- "^MT-"

# Defaults for any run with no matching library in qc_thresholds.csv
DEFAULT_QC <- list(nFeature_min = 200, nFeature_max = 5000,
                   mt_max = 15, mt_min = NA)

# -----------------------------------------------------------------------------
# 4. Normalization / integration / clustering / UMAP
# -----------------------------------------------------------------------------
N_VARIABLE_FEATURES <- 2000
N_PCS_COMPUTE       <- 50
DIMS_USE            <- 1:30          # Seurat default used for anchors, neighbors, UMAP
CLUSTER_RESOLUTIONS <- seq(0.2, 2.0, by = 0.2)
UMAP_N_NEIGHBORS    <- 40            # from the paper
UMAP_MIN_DIST       <- 0.3           # from the paper

# CCA integration: default k.weight is 100. If the smallest run has fewer
# cells than that, the script lowers k.weight automatically and logs it.
CCA_K_WEIGHT_DEFAULT <- 100

# Optional: SRR IDs to use as integration references (NULL = all pairs,
# which is what Seurat v3 did by default). Setting e.g. the MNC runs as
# references is much faster but is a deviation from the paper.
INTEGRATION_REFERENCE_RUNS <- NULL

# -----------------------------------------------------------------------------
# 5. BoneMarrowMap
# -----------------------------------------------------------------------------
BMM_BATCH_VAR      <- "Run"      # batch variable corrected during mapping
BMM_MAD_THRESHOLD  <- 2.5        # README default
BMM_BROAD_COL      <- "predicted_CellType_Broad"

# -----------------------------------------------------------------------------
# 6. Differential expression
# -----------------------------------------------------------------------------
SEED                      <- 8     # same seed as the original DE script
DE_DEFAULT_CELLS_PER_DONOR <- 150  # per cell-type comparison (override per type)
DE_GLOBAL_CELLS_PER_DONOR  <- 3000 # global FPD vs Healthy
DE_MIN_CELLS_PER_GROUP     <- 100  # after downsampling, both groups need >= this
DE_TOP_GENES_TO_LABEL      <- 30
DE_PADJ_CUTOFF             <- 0.05 # only used for the summary count

# -----------------------------------------------------------------------------
# 7. Runtime settings
# -----------------------------------------------------------------------------
# Large objects cross the 'future' globals limit inside IntegrateLayers.
options(future.globals.maxSize = 100 * 1024^3)   # 100 GB
# Keep everything sequential: parallel 'future' workers each copy the whole
# object and would multiply memory use.
options(stringsAsFactors = FALSE, warn = 1)

###############################################################################
# Helper functions
###############################################################################

# Timestamped log line (stdout is redirected to Output/logs by the run_*.sh)
log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush.console()
}

# Log R memory use. 'max used' is the peak since the last gc(reset = TRUE).
log_mem <- function(label = "") {
  g <- gc()
  log_msg("MEMORY %s | current R heap: %.1f GB | peak since reset: %.1f GB",
          label, sum(g[, 2]) / 1024, sum(g[, 6]) / 1024)
}

# Log disk usage of the Output folder
log_disk <- function() {
  du <- tryCatch(system2("du", c("-sh", shQuote(OUT_DIR)), stdout = TRUE),
                 error = function(e) "du not available")
  log_msg("DISK Output folder: %s", paste(du, collapse = " "))
}

# Create folders if needed
ensure_dirs <- function(...) {
  for (d in c(...)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Save a ggplot/patchwork object as PNG (primary, for viewing) and,
# optionally, PDF (for figures). Works on a headless node: PNG uses ragg if
# installed, otherwise the cairo png device; PDF uses the base pdf device.
save_plot <- function(p, dir, name, width = 7, height = 6, dpi = 150, pdf = TRUE) {
  ensure_dirs(dir)
  png_file <- file.path(dir, paste0(name, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_file, p, width = width, height = height, dpi = dpi,
                    device = ragg::agg_png, bg = "white", limitsize = FALSE)
  } else {
    ggplot2::ggsave(png_file, p, width = width, height = height, dpi = dpi,
                    type = "cairo", bg = "white", limitsize = FALSE)
  }
  if (pdf) {
    ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = width,
                    height = height, device = grDevices::pdf, limitsize = FALSE)
  }
  invisible(png_file)
}

# Filesystem-safe version of a label (cell type names contain spaces, '/', '+')
safe_name <- function(x) {
  x <- gsub("\\+", "pos", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# Write sessionInfo for reproducibility
write_session_info <- function(script_name) {
  ensure_dirs(LOG_DIR)
  writeLines(capture.output(sessionInfo()),
             file.path(LOG_DIR, paste0(script_name, "_sessionInfo.txt")))
}

# Delete an intermediate checkpoint (only when the setting above allows it)
delete_checkpoint <- function(path) {
  if (DELETE_INTERMEDIATE_CHECKPOINTS && file.exists(path)) {
    unlink(path)
    log_msg("Deleted intermediate checkpoint %s", basename(path))
  }
}

# Adjusted Rand index between two labelings (no extra packages needed)
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  comb2 <- function(v) sum(v * (v - 1) / 2)
  n        <- sum(tab)
  sum_ij   <- comb2(tab)
  sum_i    <- comb2(rowSums(tab))
  sum_j    <- comb2(colSums(tab))
  expected <- sum_i * sum_j / (n * (n - 1) / 2)
  (sum_ij - expected) / ((sum_i + sum_j) / 2 - expected)
}

# Normalized mutual information (geometric-mean normalization)
normalized_mutual_info <- function(x, y) {
  p  <- table(x, y) / length(x)
  px <- rowSums(p); py <- colSums(p)
  nz <- p > 0
  mi <- sum(p[nz] * log(p[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

# A long, fixed categorical palette for cell types / donors.
# Colors are assigned in a fixed order to the sorted labels so the same
# label always gets the same color across every plot in the pipeline.
fixed_palette <- function(labels) {
  base <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
            "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78",
            "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7",
            "#DBDB8D", "#9EDAE5", "#393B79", "#637939", "#8C6D31", "#843C39",
            "#7B4173", "#3182BD", "#E6550D", "#31A354", "#756BB1", "#636363")
  labels <- sort(unique(as.character(labels)))
  cols <- if (length(labels) <= length(base)) base[seq_along(labels)] else
    grDevices::colorRampPalette(base)(length(labels))
  stats::setNames(cols, labels)
}

CONDITION_COLORS <- c(FPD = "#C0392B", Healthy = "#2471A3")

# Natural donor order for plots: FPD1, FPD2, ..., FPD12, HD19002, ...
donor_levels <- function(d) {
  d <- unique(as.character(d))
  d[order(sub("[0-9]+$", "", d), suppressWarnings(as.numeric(sub("^[^0-9]*", "", d))))]
}

# "RNA_snn_res.1" -> "1.0" (so files sort and read consistently)
res_label <- function(col) sprintf("%.1f", as.numeric(sub("^RNA_snn_res\\.", "", col)))
