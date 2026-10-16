###############################################################################
# 04_differential_expression.R
#
# What this script does
#   FPD vs Healthy differential expression, following the paper:
#     FindMarkers(default Wilcoxon test, logfc.threshold = -Inf, min.pct = -Inf,
#                 min.diff.pct = -Inf)  -> every gene is tested and reported
#   To control for donor effects, each comparison is first downsampled to at
#   most N cells per donor (fixed seed). A comparison runs only if both groups
#   still have >= 100 cells.
#
#   Comparisons
#     - Global: all cells, FPD vs Healthy, N = 3000 per donor. Cells failing
#       BoneMarrowMap mapping QC ARE included here.
#     - One per BoneMarrowMap Broad cell type: <type>_FPD vs <type>_Healthy,
#       N from metadata/downsample_overrides.csv (default 150). Cells failing
#       mapping QC have no Broad label and are NOT included.
#
#   Per comparison (folder Output/04_de/<comparison>/)
#     - DE table (.tsv), positive avg_log2FC = higher in FPD
#     - Average-expression scatter, top 30 genes labeled
#     - Cells per donor after downsampling
#     - Cell proportion plot (global: composition; cell type: abundance/donor)
#   Summary: Output/04_de/de_summary.csv and de_skipped.csv
#
# Inputs : Output/checkpoints/03_FPD_BM_annotated_FINAL.rds
#          metadata/downsample_overrides.csv  (copy + edit the TEMPLATE from 03)
# Run    : Rscript 04_differential_expression.R
###############################################################################

get_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
source(file.path(get_script_dir(), "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

ensure_dirs(DE_DIR, LOG_DIR)
log_msg("=== 04_differential_expression.R started ===")
gc(reset = TRUE)

if (!requireNamespace("presto", quietly = TRUE)) {
  log_msg("WARNING: 'presto' is not installed; FindMarkers will use the much slower Wilcoxon implementation")
}

###############################################################################
# Step 0. Inputs
###############################################################################
if (!file.exists(FINAL_OBJ)) stop("Final object not found: ", FINAL_OBJ, " (run script 03)")
if (!file.exists(DOWNSAMPLE_OVERRIDES)) {
  stop("\n", DOWNSAMPLE_OVERRIDES, " not found.\n",
       "Copy ", DOWNSAMPLE_TEMPLATE, " to that name, edit the CellsPerDonor column\n",
       "(e.g. 30 or 50 for rare types, as in the original analysis), then rerun.")
}

obj <- readRDS(FINAL_OBJ)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj, search = "data")) != 1) obj <- JoinLayers(obj)
log_msg("Loaded %d genes x %d cells", nrow(obj), ncol(obj))

md <- obj@meta.data
broad <- md[[BMM_BROAD_COL]]
broad_levels <- sort(unique(na.omit(broad)))
broad_cols <- fixed_palette(broad_levels)

overrides <- read.csv(DOWNSAMPLE_OVERRIDES, stringsAsFactors = FALSE)
if (!all(c("CellType_Broad", "CellsPerDonor") %in% colnames(overrides))) {
  stop("downsample_overrides.csv needs columns CellType_Broad and CellsPerDonor")
}
overrides$CellsPerDonor <- as.integer(overrides$CellsPerDonor)
unknown <- setdiff(overrides$CellType_Broad, broad_levels)
if (length(unknown) > 0) stop("Unknown CellType_Broad in overrides: ", paste(unknown, collapse = ", "))
if (any(is.na(overrides$CellsPerDonor) | overrides$CellsPerDonor < 1)) {
  stop("CellsPerDonor must be a positive whole number for every row")
}
n_per_donor_for <- function(type) {
  v <- overrides$CellsPerDonor[overrides$CellType_Broad == type]
  if (length(v) == 1) v else DE_DEFAULT_CELLS_PER_DONOR
}

# Normalized data (log1p CP10K) used for the average-expression scatter
norm_data <- LayerData(obj, assay = "RNA", layer = "data")

###############################################################################
# Helper functions
###############################################################################

# Downsample to at most n cells per donor. Donors with fewer cells keep all.
# Donors are processed in sorted order with one fixed seed -> reproducible.
downsample_per_donor <- function(cells, donors, n, seed = SEED) {
  set.seed(seed)
  keep <- character(0)
  for (d in sort(unique(donors))) {
    dc <- cells[donors == d]
    keep <- c(keep, if (length(dc) > n) sample(dc, n) else dc)
  }
  keep
}

# log1p(mean(expm1(x))) per gene = what Seurat's AverageExpression reports
avg_expression <- function(cells) {
  x <- norm_data[, cells, drop = FALSE]
  x@x <- expm1(x@x)
  log1p(Matrix::rowMeans(x))
}

run_comparison <- function(name, candidate_cells, n_per_donor, ident1, ident2, type = NULL) {
  out_dir <- file.path(DE_DIR, safe_name(name))

  donors <- md[candidate_cells, "Donor"]
  sel <- downsample_per_donor(candidate_cells, donors, n_per_donor)
  cond <- md[sel, "Condition"]
  cells1 <- sel[cond == "FPD"]
  cells2 <- sel[cond == "Healthy"]

  info <- data.frame(Comparison = name, ident.1 = ident1, ident.2 = ident2,
                     CellsPerDonor = n_per_donor,
                     candidate_FPD = sum(md[candidate_cells, "Condition"] == "FPD"),
                     candidate_Healthy = sum(md[candidate_cells, "Condition"] == "Healthy"),
                     n_FPD = length(cells1), n_Healthy = length(cells2),
                     FPD_donors = length(unique(md[cells1, "Donor"])),
                     Healthy_donors = length(unique(md[cells2, "Donor"])))

  if (length(cells1) < DE_MIN_CELLS_PER_GROUP || length(cells2) < DE_MIN_CELLS_PER_GROUP) {
    info$reason <- sprintf("fewer than %d cells in a group after downsampling", DE_MIN_CELLS_PER_GROUP)
    log_msg("SKIP %-35s FPD=%d Healthy=%d", name, length(cells1), length(cells2))
    return(list(skipped = info))
  }
  ensure_dirs(out_dir)
  log_msg("RUN  %-35s FPD=%d (%d donors) Healthy=%d (%d donors), N/donor=%d", name,
          length(cells1), info$FPD_donors, length(cells2), info$Healthy_donors, n_per_donor)

  # ---- Differential expression (all genes) ----
  de <- FindMarkers(obj, ident.1 = cells1, ident.2 = cells2, assay = "RNA",
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf,
                    verbose = FALSE)
  de <- de[order(de$p_val_adj, de$p_val, -abs(de$avg_log2FC)), ]
  avg1 <- avg_expression(cells1); avg2 <- avg_expression(cells2)
  de_out <- data.frame(gene = rownames(de), de, row.names = NULL, check.names = FALSE)
  de_out$avg_expr_FPD     <- avg1[de_out$gene]
  de_out$avg_expr_Healthy <- avg2[de_out$gene]
  file_stub <- paste("DE", ident1, ident2, sep = ".")
  write.table(de_out, file.path(out_dir, paste0(safe_name(file_stub), ".Markers.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- Average-expression scatter ----
  sc <- data.frame(gene = names(avg1), FPD = avg1, Healthy = avg2[names(avg1)])
  sc$significant <- sc$gene %in% de_out$gene[de_out$p_val_adj < DE_PADJ_CUTOFF]
  top <- head(de_out$gene, DE_TOP_GENES_TO_LABEL)
  p_sc <- ggplot(sc, aes(FPD, Healthy)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
    geom_point(data = sc[!sc$significant, ], colour = "grey70", size = 0.6) +
    geom_point(data = sc[sc$significant, ], colour = "#C0392B", size = 0.8) +
    geom_text_repel(data = sc[sc$gene %in% top, ], aes(label = gene), size = 2.6,
                    max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey40") +
    labs(x = sprintf("%s  [log1p mean expression]", ident1),
         y = sprintf("%s  [log1p mean expression]", ident2),
         title = sprintf("DE of %s vs %s", ident1, ident2),
         subtitle = sprintf("%d vs %d cells; red = adj. p < %.2f; top %d genes labeled",
                            length(cells1), length(cells2), DE_PADJ_CUTOFF, DE_TOP_GENES_TO_LABEL)) +
    theme_classic()
  save_plot(p_sc, out_dir, paste0(safe_name(file_stub), ".avgExpScatter"), width = 6.5, height = 6.5)

  # ---- Cells per donor after downsampling ----
  dc <- md[c(cells1, cells2), ] %>% count(Donor, Condition)
  p_dc <- ggplot(dc, aes(reorder(Donor, n), n, fill = Condition)) +
    geom_col(width = 0.8) + coord_flip() +
    geom_hline(yintercept = n_per_donor, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = CONDITION_COLORS) +
    labs(x = NULL, y = "Cells used", title = sprintf("%s: cells per donor", name),
         subtitle = sprintf("Dashed line = cap of %d cells per donor", n_per_donor)) +
    theme_classic()
  save_plot(p_dc, out_dir, paste0(safe_name(file_stub), ".DonorCounts"), width = 6.5, height = 5)
  write.csv(dc, file.path(out_dir, paste0(safe_name(file_stub), ".DonorCounts.csv")), row.names = FALSE)

  # ---- Cell proportion plot ----
  if (is.null(type)) {
    # Global: Broad composition of the cells used, per condition and donor
    comp <- md[c(cells1, cells2), ] %>%
      mutate(Broad = ifelse(is.na(.data[[BMM_BROAD_COL]]), "Mapping QC fail", .data[[BMM_BROAD_COL]])) %>%
      count(Donor, Condition, Broad) %>% group_by(Donor) %>% mutate(frac = n / sum(n)) %>% ungroup()
    comp$Broad <- factor(comp$Broad, levels = c(broad_levels, "Mapping QC fail"))
    comp$Donor <- factor(comp$Donor, levels = donor_levels(comp$Donor))
    p_prop <- ggplot(comp, aes(Donor, frac, fill = Broad)) +
      geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
      facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = c(broad_cols, "Mapping QC fail" = "#D9D9D9")) +
      scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
      labs(x = NULL, y = "% of cells used", fill = "Broad cell type",
           title = "Cell-type composition of the cells in this comparison") +
      theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p_prop, out_dir, paste0(safe_name(file_stub), ".cellProportions"), width = 11, height = 6)
  } else {
    # Cell type: abundance of this type in each donor (all QC-pass cells, before downsampling)
    ab <- md[!is.na(md[[BMM_BROAD_COL]]), ] %>%
      group_by(Donor, Condition) %>%
      summarise(pct = 100 * mean(.data[[BMM_BROAD_COL]] == type), .groups = "drop")
    set.seed(SEED)
    ab$x <- as.numeric(factor(ab$Condition, levels = c("FPD", "Healthy"))) +
      stats::runif(nrow(ab), -0.12, 0.12)
    med <- ab %>% group_by(Condition) %>% summarise(m = median(pct), .groups = "drop") %>%
      mutate(x = as.numeric(factor(Condition, levels = c("FPD", "Healthy"))))
    p_prop <- ggplot(ab, aes(x, pct, colour = Condition)) +
      geom_segment(data = med, aes(x = x - 0.25, xend = x + 0.25, y = m, yend = m),
                   colour = "grey30", linewidth = 0.8) +
      geom_point(size = 2.5) +
      geom_text_repel(aes(label = Donor), size = 2.5, colour = "grey30") +
      scale_x_continuous(breaks = 1:2, labels = c("FPD", "Healthy"), limits = c(0.5, 2.5)) +
      scale_colour_manual(values = CONDITION_COLORS) +
      labs(x = NULL, y = sprintf("%% of donor's QC-pass cells that are %s", type),
           title = sprintf("Abundance of %s per donor", type),
           subtitle = "Libraries pooled per donor; bar = median") +
      theme_classic() + theme(legend.position = "none")
    save_plot(p_prop, out_dir, paste0(safe_name(file_stub), ".cellProportions"), width = 5.5, height = 5.5)
  }

  info$n_genes_tested   <- nrow(de_out)
  info$n_sig_padj       <- sum(de_out$p_val_adj < DE_PADJ_CUTOFF)
  info$n_sig_up_in_FPD  <- sum(de_out$p_val_adj < DE_PADJ_CUTOFF & de_out$avg_log2FC > 0)
  info$n_sig_down_in_FPD <- sum(de_out$p_val_adj < DE_PADJ_CUTOFF & de_out$avg_log2FC < 0)
  info$output_folder    <- out_dir
  list(done = info)
}

###############################################################################
# Step 1. Run all comparisons
###############################################################################
results <- list()

# Global FPD vs Healthy (all cells, including mapping-QC-fail cells)
results[["Global"]] <- run_comparison(
  name = "Global_FPD_vs_Healthy", candidate_cells = rownames(md),
  n_per_donor = DE_GLOBAL_CELLS_PER_DONOR, ident1 = "FPD", ident2 = "Healthy")
invisible(gc())

# One comparison per Broad cell type
for (type in broad_levels) {
  cand <- rownames(md)[!is.na(broad) & broad == type]
  results[[type]] <- run_comparison(
    name = paste0(type, "_FPD_vs_Healthy"), candidate_cells = cand,
    n_per_donor = n_per_donor_for(type),
    ident1 = paste0(type, "_FPD"), ident2 = paste0(type, "_Healthy"), type = type)
  invisible(gc())
}

###############################################################################
# Step 2. Summary tables
###############################################################################
done    <- bind_rows(lapply(results, `[[`, "done"))
skipped <- bind_rows(lapply(results, `[[`, "skipped"))
write.csv(done, file.path(DE_DIR, "de_summary.csv"), row.names = FALSE)
if (nrow(skipped) > 0) {
  write.csv(skipped, file.path(DE_DIR, "de_skipped.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(note = "no comparisons were skipped"), file.path(DE_DIR, "de_skipped.csv"),
            row.names = FALSE)
}
log_msg("Comparisons run: %d | skipped: %d", nrow(done), nrow(skipped))
if (nrow(done) > 0) print(done[, c("Comparison", "CellsPerDonor", "n_FPD", "n_Healthy", "n_sig_padj")])

# Overview: number of significant genes per comparison
if (nrow(done) > 0) {
  ov <- done %>% select(Comparison, n_sig_up_in_FPD, n_sig_down_in_FPD) %>%
    tidyr::pivot_longer(-Comparison, names_to = "Direction", values_to = "n") %>%
    mutate(Direction = ifelse(Direction == "n_sig_up_in_FPD", "Higher in FPD", "Lower in FPD"))
  p_ov <- ggplot(ov, aes(reorder(Comparison, n), n, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) + coord_flip() +
    scale_fill_manual(values = c("Higher in FPD" = "#C0392B", "Lower in FPD" = "#2471A3")) +
    labs(x = NULL, y = sprintf("Genes with adj. p < %.2f", DE_PADJ_CUTOFF),
         title = "Significant genes per comparison") +
    theme_classic()
  save_plot(p_ov, DE_DIR, "de_overview_significant_genes", width = 9,
            height = max(4, 0.35 * nrow(done) + 2))
}

log_mem("end of 04")
log_disk()
write_session_info("04_differential_expression")
log_msg("=== 04_differential_expression.R finished ===")
