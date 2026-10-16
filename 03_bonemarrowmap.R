###############################################################################
# 03_bonemarrowmap.R
#
# What this script does
#   1. Maps every cell onto the BoneMarrowMap reference atlas, following the
#      README Quick Start (https://github.com/andygxzeng/BoneMarrowMap):
#        map_Query -> calculate_MappingError -> predict_CellTypes ->
#        predict_Pseudotime -> plot_Projection_byDonor -> save_ProjectionResults
#   2. Copies the results into the integrated object's metadata:
#        mapping_error_score, mapping_error_QC (Pass/Fail),
#        initial_CellType (all cells), predicted_CellType (NA if QC Fail),
#        predicted_CellType_Broad, predicted_CellType_No_Abbreviations,
#        initial_Pseudotime, predicted_Pseudotime, and the projected UMAP.
#   3. Compares each clustering resolution with the Broad labels
#      (row-normalized heatmap + adjusted Rand index / normalized mutual info)
#      so you can pick the resolution that best matches BoneMarrowMap.
#   4. Plots Broad cell-type proportions by donor, condition and library type.
#   5. Writes the cell counts needed to choose DE downsampling numbers, and a
#      downsample_overrides_TEMPLATE.csv for you to edit before script 04.
#
# Input  : Output/checkpoints/02_integrated_clustered.rds
#          BoneMarrowMap_ref/BoneMarrowMap_SymphonyReference.rds
#          BoneMarrowMap_ref/BoneMarrowMap_uwot_model.uwot
# Outputs: Output/03_bonemarrowmap/, Output/checkpoints/03_FPD_BM_annotated_FINAL.rds
#
# Run    : Rscript 03_bonemarrowmap.R
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
  library(symphony)
  library(BoneMarrowMap)
})

DONOR_PROJ_DIR <- file.path(BMM_OUT, "projection_by_donor")
RES_DIR        <- file.path(BMM_OUT, "clusters_vs_BoneMarrowMap")
PROP_DIR       <- file.path(BMM_OUT, "celltype_proportions")
ensure_dirs(BMM_OUT, DONOR_PROJ_DIR, RES_DIR, PROP_DIR, LOG_DIR, META_DIR)
log_msg("=== 03_bonemarrowmap.R started ===")
gc(reset = TRUE)
set.seed(SEED)

###############################################################################
# Step 0. Load the integrated object and the reference
###############################################################################
if (file.exists(CKPT_02)) {
  obj <- readRDS(CKPT_02)
} else if (file.exists(FINAL_OBJ)) {
  log_msg("02 checkpoint not found (already deleted); re-running on the FINAL object")
  obj <- readRDS(FINAL_OBJ)
} else {
  stop("Neither ", CKPT_02, " nor ", FINAL_OBJ, " exists. Run script 02 first.")
}
log_msg("Loaded %d genes x %d cells", nrow(obj), ncol(obj))

for (f in c(BMM_REF_RDS, BMM_UWOT)) if (!file.exists(f)) stop("Missing reference file: ", f)
ref <- readRDS(BMM_REF_RDS)
ref$save_uwot_path <- BMM_UWOT          # README: set uwot path for UMAP projection
log_msg("Loaded BoneMarrowMap reference")

# Remove any earlier BoneMarrowMap columns (makes this script safe to rerun)
bmm_cols_old <- grep("^(mapping_error|initial_|predicted_)", colnames(obj@meta.data), value = TRUE)
if (length(bmm_cols_old) > 0) {
  obj@meta.data[bmm_cols_old] <- NULL
  log_msg("Removed old BoneMarrowMap columns from a previous run: %s", paste(bmm_cols_old, collapse = ", "))
}

###############################################################################
# Step 1. Map query cells (README Quick Start)
###############################################################################
# A lightweight query object: raw counts + metadata only. map_Query reads the
# RNA 'counts' layer and normalizes it itself (log CP10K), as the README expects.
query <- CreateSeuratObject(
  counts    = LayerData(obj, assay = "RNA", layer = "counts"),
  meta.data = obj@meta.data[, c("Run", "Donor", "Condition", "Sex", "LibraryType",
                                "nCount_RNA", "nFeature_RNA", "percent.mt")])

log_msg("Running map_Query (batch variable = %s)", BMM_BATCH_VAR)
query <- map_Query(query = query, ref_obj = ref, vars = BMM_BATCH_VAR)

query <- calculate_MappingError(query, reference = ref, MAD_threshold = BMM_MAD_THRESHOLD)
log_msg("Mapping QC: %d Pass, %d Fail", sum(query$mapping_error_QC == "Pass"),
        sum(query$mapping_error_QC == "Fail"))

mapqc_plots <- plot_MappingErrorQC(query)
save_plot(wrap_plots(mapqc_plots, ncol = 2), BMM_OUT, "mapping_error_QC", width = 12,
          height = 10, pdf = FALSE)

query <- predict_CellTypes(query_obj = query, ref_obj = ref,
                           final_label = "predicted_CellType")
query <- predict_Pseudotime(query_obj = query, ref_obj = ref,
                            final_label = "predicted_Pseudotime")
log_mem("after mapping")

# Per-donor density plots on the reference UMAP. The package's own file saving
# builds a relative path (fails with absolute paths), so plots are returned
# and saved here instead.
donor_plots <- plot_Projection_byDonor(query_obj = query, batch_key = "Donor",
                                       ref_obj = ref, saveplot = FALSE)
for (d in names(donor_plots)) {
  save_plot(donor_plots[[d]], DONOR_PROJ_DIR, paste0("density_", safe_name(d), "_projectedUMAP"),
            width = 6, height = 4, pdf = TRUE)
}
save_plot(wrap_plots(donor_plots, ncol = 5), DONOR_PROJ_DIR, "ALL_donors_projectedUMAP",
          width = 22, height = 12, pdf = FALSE)

save_ProjectionResults(query_obj = query,
                       file_name = file.path(BMM_OUT, "querydata_projected_labeled.csv"))
log_msg("Saved per-cell projection results")

###############################################################################
# Step 2. Copy results into the integrated object
###############################################################################
if (!identical(colnames(query), colnames(obj))) stop("Cell barcodes differ between query and object")
new_cols <- setdiff(colnames(query@meta.data), colnames(obj@meta.data))
new_cols <- setdiff(new_cols, c("orig.ident", "nCount_RNA", "nFeature_RNA"))
obj <- AddMetaData(obj, metadata = query@meta.data[colnames(obj), new_cols, drop = FALSE])
obj[["umap_projected"]] <- CreateDimReducObject(
  embeddings = Embeddings(query, reduction = "umap_projected")[colnames(obj), ],
  key = "umapprojected_", assay = "RNA")
log_msg("Added metadata columns: %s", paste(new_cols, collapse = ", "))

ref_umap <- data.frame(ref$umap$embedding)
colnames(ref_umap) <- c("x", "y")
proj_umap <- data.frame(Embeddings(query, "umap_projected"))
rm(query); invisible(gc())

if (!BMM_BROAD_COL %in% colnames(obj@meta.data)) {
  stop(BMM_BROAD_COL, " was not created. Is the reference object up to date?")
}
broad_levels <- sort(unique(na.omit(obj[[BMM_BROAD_COL, drop = TRUE]])))
broad_cols   <- fixed_palette(broad_levels)
log_msg("Broad cell types (%d): %s", length(broad_levels), paste(broad_levels, collapse = " | "))

# Plotting-only column: QC-failing cells shown as their own grey group
FAIL_LABEL <- "Mapping QC fail"
obj$Broad_for_plot <- ifelse(is.na(obj[[BMM_BROAD_COL, drop = TRUE]]), FAIL_LABEL,
                             obj[[BMM_BROAD_COL, drop = TRUE]])
obj$Broad_for_plot <- factor(obj$Broad_for_plot, levels = c(broad_levels, FAIL_LABEL))
broad_cols_fail <- c(broad_cols, setNames("#D9D9D9", FAIL_LABEL))

###############################################################################
# Step 3. Annotation plots
###############################################################################
p_broad <- DimPlot(obj, reduction = "umap", group.by = "Broad_for_plot", raster = TRUE,
                   label = TRUE, repel = TRUE, label.size = 3, cols = broad_cols_fail) +
  ggtitle("BoneMarrowMap Broad cell types on integrated UMAP")
save_plot(p_broad, BMM_OUT, "umap_BoneMarrowMap_Broad", width = 11, height = 8)

p_fine <- DimPlot(obj, reduction = "umap", group.by = "predicted_CellType", raster = TRUE,
                  label = TRUE, repel = TRUE, label.size = 2.2) + NoLegend() +
  ggtitle("BoneMarrowMap fine cell types (QC-pass cells; failing cells grey)")
save_plot(p_fine, BMM_OUT, "umap_BoneMarrowMap_fine", width = 11, height = 10)

p_qc <- DimPlot(obj, reduction = "umap", group.by = "mapping_error_QC", raster = TRUE,
                cols = c(Fail = "#C0392B", Pass = "#BDBDBD"), order = "Fail") +
  ggtitle("Mapping error QC on integrated UMAP")
save_plot(p_qc, BMM_OUT, "umap_mapping_error_QC", width = 8, height = 6.5)

p_pt <- FeaturePlot(obj, features = "predicted_Pseudotime", reduction = "umap", raster = TRUE) +
  ggtitle("BoneMarrowMap predicted pseudotime")
save_plot(p_pt, BMM_OUT, "umap_predicted_Pseudotime", width = 8, height = 6.5)

# Query cells on the reference UMAP, colored by Broad label
proj_df <- data.frame(x = proj_umap[, 1], y = proj_umap[, 2],
                      Broad = obj$Broad_for_plot, Condition = obj$Condition)
set.seed(SEED)
bg <- ref_umap[sample(nrow(ref_umap), min(60000, nrow(ref_umap))), ]
p_proj <- ggplot(proj_df[sample(nrow(proj_df)), ], aes(x, y)) +
  geom_point(data = bg, colour = "#EBEBEB", size = 0.05) +
  geom_point(aes(colour = Broad), size = 0.05, alpha = 0.6) +
  scale_colour_manual(values = broad_cols_fail) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~Condition) + theme_void() +
  ggtitle("Cells projected onto the BoneMarrowMap reference UMAP (grey = reference)")
save_plot(p_proj, BMM_OUT, "reference_UMAP_projection_by_condition", width = 14, height = 6.5, pdf = FALSE)

# Mapping QC pass rate by run and by library type
qc_rate <- obj@meta.data %>%
  group_by(Run, Donor, Condition, LibraryType) %>%
  summarise(n_cells = n(), n_fail = sum(mapping_error_QC == "Fail"),
            pct_fail = round(100 * n_fail / n_cells, 1), .groups = "drop")
write.csv(qc_rate, file.path(BMM_OUT, "mapping_QC_fail_rate_by_run.csv"), row.names = FALSE)
p_rate <- ggplot(qc_rate, aes(reorder(Run, pct_fail), pct_fail, fill = LibraryType)) +
  geom_col(width = 0.8) + coord_flip() +
  labs(x = NULL, y = "% cells failing mapping QC", title = "Mapping QC failure by run") +
  theme_classic()
save_plot(p_rate, BMM_OUT, "mapping_QC_fail_rate_by_run", width = 8, height = 8)

###############################################################################
# Step 4. Clusters vs BoneMarrowMap Broad labels, per resolution
###############################################################################
res_cols <- grep("^RNA_snn_res\\.", colnames(obj@meta.data), value = TRUE)
res_cols <- res_cols[order(as.numeric(sub("RNA_snn_res\\.", "", res_cols)))]
pass <- !is.na(obj[[BMM_BROAD_COL, drop = TRUE]])
broad_pass <- obj[[BMM_BROAD_COL, drop = TRUE]][pass]

agreement <- list()
for (rc in res_cols) {
  res <- res_label(rc)
  cl  <- obj[[rc, drop = TRUE]]

  agreement[[rc]] <- data.frame(
    Resolution = as.numeric(res), n_clusters = length(unique(cl)),
    ARI = adjusted_rand_index(cl[pass], broad_pass),
    NMI = normalized_mutual_info(cl[pass], broad_pass))

  # Row-normalized contingency: what fraction of each cluster is each label
  tab <- as.data.frame(table(Cluster = cl[pass], Broad = broad_pass)) %>%
    group_by(Cluster) %>% mutate(n_cluster = sum(Freq), frac = Freq / pmax(n_cluster, 1)) %>%
    ungroup()
  majority <- tab %>% group_by(Cluster) %>% slice_max(frac, n = 1, with_ties = FALSE) %>%
    transmute(Cluster, n_cells_QCpass = n_cluster, majority_Broad = Broad,
              majority_fraction = round(frac, 3))
  write.csv(majority, file.path(RES_DIR, sprintf("cluster_majority_label_res_%s.csv", res)),
            row.names = FALSE)

  p_heat <- ggplot(tab, aes(Broad, Cluster, fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(frac >= 0.10, round(100 * frac), "")), size = 2.3) +
    scale_fill_gradient(low = "#FFFFFF", high = "#08306B", limits = c(0, 1),
                        name = "Fraction of\ncluster") +
    labs(x = "BoneMarrowMap Broad label", y = "Cluster",
         title = sprintf("Resolution %s: cluster composition (numbers = %%, shown if >= 10%%)", res)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

  p_cl <- DimPlot(obj, reduction = "umap", group.by = rc, label = TRUE, label.size = 3,
                  repel = TRUE, raster = TRUE) + NoLegend() +
    ggtitle(sprintf("Clusters, resolution %s", res))
  p_bd <- DimPlot(obj, reduction = "umap", group.by = "Broad_for_plot", raster = TRUE,
                  cols = broad_cols_fail) + ggtitle("BoneMarrowMap Broad") +
    theme(legend.text = element_text(size = 7))
  composite <- (p_cl | p_bd) / p_heat + plot_layout(heights = c(1, 1.1))
  save_plot(composite, RES_DIR, sprintf("res_%s_clusters_vs_Broad", res),
            width = 17, height = 15, pdf = FALSE)
}
agreement <- bind_rows(agreement)
write.csv(agreement, file.path(RES_DIR, "agreement_ARI_NMI_by_resolution.csv"), row.names = FALSE)
print(agreement)

agree_long <- agreement %>% pivot_longer(c(ARI, NMI), names_to = "Metric", values_to = "Value")
p_agree <- ggplot(agree_long, aes(Resolution, Value, colour = Metric)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  geom_text(data = agreement, aes(Resolution, 1.02, label = n_clusters), inherit.aes = FALSE,
            size = 3, colour = "grey30") +
  scale_colour_manual(values = c(ARI = "#1F77B4", NMI = "#FF7F0E")) +
  scale_x_continuous(breaks = agreement$Resolution) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(y = "Agreement with BoneMarrowMap Broad labels",
       title = "Cluster vs BoneMarrowMap agreement (numbers = clusters)",
       subtitle = "Mapping-QC-pass cells only; higher = better match") +
  theme_classic()
save_plot(p_agree, RES_DIR, "agreement_ARI_NMI_by_resolution", width = 8, height = 5)

###############################################################################
# Step 5. Broad cell-type proportions (overall)
###############################################################################
md_pass <- obj@meta.data[pass, ]
md_pass$Broad <- factor(md_pass[[BMM_BROAD_COL]], levels = broad_levels)
md_pass$Donor <- factor(md_pass$Donor, levels = donor_levels(md_pass$Donor))

prop_donor <- md_pass %>% count(Donor, Condition, Broad) %>%
  group_by(Donor) %>% mutate(frac = n / sum(n)) %>% ungroup()
p_prop_donor <- ggplot(prop_donor, aes(Donor, frac, fill = Broad)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = broad_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "% of donor's QC-pass cells", fill = "Broad cell type",
       title = "BoneMarrowMap Broad composition per donor (all libraries pooled)") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_prop_donor, PROP_DIR, "proportions_by_donor", width = 12, height = 6)

prop_lib <- md_pass %>% count(Donor, Condition, LibraryType, Broad) %>%
  group_by(Donor, LibraryType) %>% mutate(frac = n / sum(n)) %>% ungroup()
p_prop_lib <- ggplot(prop_lib, aes(Donor, frac, fill = Broad)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  facet_grid(LibraryType ~ Condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = broad_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "% of cells in library type", fill = "Broad cell type",
       title = "Composition split by library type (enrichment changes composition!)") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_prop_lib, PROP_DIR, "proportions_by_donor_and_librarytype", width = 13, height = 10)

prop_cond <- md_pass %>% count(Condition, LibraryType, Broad) %>%
  group_by(Condition, LibraryType) %>% mutate(frac = n / sum(n)) %>% ungroup()
p_prop_cond <- ggplot(prop_cond, aes(Condition, frac, fill = Broad)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  facet_wrap(~LibraryType, nrow = 1) +
  scale_fill_manual(values = broad_cols) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "% of cells", fill = "Broad cell type",
       title = "Composition by condition, within each library type") +
  theme_classic()
save_plot(p_prop_cond, PROP_DIR, "proportions_by_condition_and_librarytype", width = 11, height = 6)
write.csv(prop_lib, file.path(PROP_DIR, "proportions_by_donor_and_librarytype.csv"), row.names = FALSE)

###############################################################################
# Step 6. Cell counts for DE planning + downsampling template
###############################################################################
counts_long <- md_pass %>% count(Broad, Donor, Condition, name = "n_cells")
write.csv(counts_long, file.path(BMM_OUT, "Broad_celltype_counts_by_donor.csv"), row.names = FALSE)

counts_wide <- counts_long %>%
  pivot_wider(names_from = Donor, values_from = n_cells, values_fill = 0, id_cols = Broad)
write.csv(counts_wide, file.path(BMM_OUT, "Broad_celltype_counts_by_donor_WIDE.csv"), row.names = FALSE)

# Projected cells per group if downsampled to N per donor
proj_after <- function(N) {
  counts_long %>% group_by(Broad, Condition) %>%
    summarise(n = sum(pmin(n_cells, N)), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = n, values_fill = 0)
}
tmpl <- counts_long %>% group_by(Broad) %>%
  summarise(total_FPD      = sum(n_cells[Condition == "FPD"]),
            total_Healthy  = sum(n_cells[Condition == "Healthy"]),
            FPD_donors     = sum(Condition == "FPD" & n_cells > 0),
            Healthy_donors = sum(Condition == "Healthy" & n_cells > 0), .groups = "drop")
for (N in c(30, 50, DE_DEFAULT_CELLS_PER_DONOR)) {
  pa <- proj_after(N)
  tmpl[[paste0("FPD_after_", N)]]     <- pa$FPD[match(tmpl$Broad, pa$Broad)]
  tmpl[[paste0("Healthy_after_", N)]] <- pa$Healthy[match(tmpl$Broad, pa$Broad)]
}
tmpl[is.na(tmpl)] <- 0
tmpl <- tmpl %>%
  mutate(passes_min_at_default =
           .data[[paste0("FPD_after_", DE_DEFAULT_CELLS_PER_DONOR)]] >= DE_MIN_CELLS_PER_GROUP &
           .data[[paste0("Healthy_after_", DE_DEFAULT_CELLS_PER_DONOR)]] >= DE_MIN_CELLS_PER_GROUP,
         CellsPerDonor = DE_DEFAULT_CELLS_PER_DONOR) %>%
  rename(CellType_Broad = Broad) %>%
  select(CellType_Broad, CellsPerDonor, everything())
write.csv(tmpl, DOWNSAMPLE_TEMPLATE, row.names = FALSE)
log_msg("Wrote %s -- edit CellsPerDonor, save as %s, then run script 04",
        DOWNSAMPLE_TEMPLATE, DOWNSAMPLE_OVERRIDES)

###############################################################################
# Save FINAL object, delete the 02 checkpoint
###############################################################################
obj$Broad_for_plot <- NULL
saveRDS(obj, FINAL_OBJ)
log_msg("Saved FINAL annotated object %s", FINAL_OBJ)
delete_checkpoint(CKPT_02)
log_mem("end of 03")
log_disk()
write_session_info("03_bonemarrowmap")
log_msg("=== 03_bonemarrowmap.R finished ===")
