###############################################################################
# 05_tf_activity.R
#
# TF and pathway activity, FPD vs Healthy, per progenitor / myeloid cell type.
#
# The statistical unit is the DONOR, not the cell. Cells are collapsed to
# pseudobulk (sum of raw counts per library x cell type), then limma-voom fits
# FPD vs Healthy with covariates, and decoupleR converts the per-gene statistics
# into TF activity (CollecTRI) and pathway activity (PROGENy).
#
# TWO CORRELATION STRUCTURES ARE FITTED AND BOTH REPORTED, because the right
# one is a judgement call, not a fact about the data:
#   block = "Donor"  - accounts for several libraries from the same person
#   block = "Family" - also accounts for relatives sharing a germline variant
#                      (more conservative; effectively fewer independent units)
# Neither is labelled primary. Compare them before deciding what to report.
#
# Inputs
#   Output/checkpoints/03_FPD_BM_annotated_FINAL.rds
#   resources/collectri_human.rds, resources/progeny_human_top500.rds
#   metadata/celltype_groups.csv   <- YOU fill this in (template written on run 1)
#
# Outputs (Output/05_tf_activity/)
#   BM_scRNA_<group>_<method>_results.tsv   per the shared spec, Section 8
#   BM_scRNA_<group>_<method>_perdonor.tsv
#   BM_scRNA_<group>_<method>_perm.tsv      label-shuffled null
#   pseudobulk_summary.csv, plots/
#
# Run: Rscript 05_tf_activity.R
###############################################################################

get_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
source(file.path(get_script_dir(), "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tidyr); library(tibble)
  library(limma); library(edgeR); library(decoupleR)
  library(ggplot2); library(patchwork)
})

TF_DIR   <- file.path(OUT_DIR, "05_tf_activity")
PLOT_DIR <- file.path(TF_DIR, "plots")
ensure_dirs(TF_DIR, PLOT_DIR, LOG_DIR, META_DIR)
log_msg("=== 05_tf_activity.R started ===")
invisible(gc(reset = TRUE))
set.seed(SEED)

###############################################################################
# Step 0. Inputs
###############################################################################
if (!file.exists(FINAL_OBJ)) stop("Annotated object not found: ", FINAL_OBJ, " (run script 03)")
for (f in c(COLLECTRI_HUMAN, PROGENY_HUMAN)) {
  if (!file.exists(f)) stop("Network file missing: ", f,
                            "\nSee DOWNLOADS_manifest.md - these must be fetched on a machine",
                            "\nwith internet and uploaded; decoupleR cannot download them here.")
}
collectri <- readRDS(COLLECTRI_HUMAN)
progeny   <- readRDS(PROGENY_HUMAN)

# decoupleR network column names have changed between versions; detect them.
net_cols <- function(net, kind) {
  src <- intersect(c("source", "tf", "pathway"), colnames(net))[1]
  tgt <- intersect(c("target", "gene"), colnames(net))[1]
  mor <- intersect(c("mor", "weight"), colnames(net))[1]
  if (is.na(src) || is.na(tgt) || is.na(mor))
    stop(kind, " network needs source/target/mor columns; found: ",
         paste(colnames(net), collapse = ", "))
  list(source = src, target = tgt, mor = mor)
}
cc <- net_cols(collectri, "CollecTRI")
pc <- net_cols(progeny, "PROGENy")
log_msg("CollecTRI: %d edges, %d TFs | PROGENy: %d edges, %d pathways",
        nrow(collectri), n_distinct(collectri[[cc$source]]),
        nrow(progeny), n_distinct(progeny[[pc$source]]))

obj <- readRDS(FINAL_OBJ)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj, search = "counts")) != 1) obj <- JoinLayers(obj)
md <- obj@meta.data
log_msg("Loaded %d genes x %d cells", nrow(obj), ncol(obj))

# Family column: donors sharing a germline variant are treated as one family.
# Donors with no group (all controls, and singleton patients) are their own.
if (!"Shared_variant_group" %in% colnames(md)) {
  log_msg("No Shared_variant_group column; every donor treated as its own family")
  md$Shared_variant_group <- md$Donor
}
md$Family <- ifelse(is.na(md$Shared_variant_group) | md$Shared_variant_group == "",
                    md$Donor, md$Shared_variant_group)

###############################################################################
# Step 1. Cell-type groups (you define these after seeing the labels)
###############################################################################
label_col <- if ("predicted_CellType" %in% colnames(md)) "predicted_CellType" else BMM_BROAD_COL
if (!file.exists(CELLTYPE_GROUPS_CSV)) {
  avail <- md %>%
    filter(!is.na(.data[[label_col]])) %>%
    group_by(CellType_fine = .data[[label_col]]) %>%
    summarise(CellType_Broad = names(sort(table(.data[[BMM_BROAD_COL]]), decreasing = TRUE))[1],
              n_cells = n(),
              n_donors_FPD = n_distinct(Donor[Condition == "FPD"]),
              n_donors_Healthy = n_distinct(Donor[Condition == "Healthy"]),
              .groups = "drop") %>%
    arrange(desc(n_cells)) %>%
    mutate(AnalysisGroup = "")
  write.csv(avail, CELLTYPE_GROUPS_TEMPLATE, row.names = FALSE, na = "")
  stop(paste0(
    "\n\ncelltype_groups.csv does not exist yet. A template listing every\n",
    "BoneMarrowMap label present in your data was written to:\n  ",
    CELLTYPE_GROUPS_TEMPLATE, "\n\n",
    "Fill in the AnalysisGroup column, e.g. HSC, MPP, CMP, GMP, LateMyeloid.\n",
    "Several fine labels may map to one group. Leave AnalysisGroup blank to\n",
    "exclude that label. Save as ", CELLTYPE_GROUPS_CSV, " and rerun.\n"))
}

groups <- read.csv(CELLTYPE_GROUPS_CSV, stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (!all(c("CellType_fine", "AnalysisGroup") %in% colnames(groups)))
  stop("celltype_groups.csv needs columns CellType_fine and AnalysisGroup")
groups <- groups %>% filter(!is.na(AnalysisGroup)) %>%
  mutate(across(c(CellType_fine, AnalysisGroup), trimws))
unknown <- setdiff(groups$CellType_fine, unique(md[[label_col]]))
if (length(unknown) > 0) stop("Labels in celltype_groups.csv not present in the data: ",
                              paste(unknown, collapse = ", "))
md$AnalysisGroup <- groups$AnalysisGroup[match(md[[label_col]], groups$CellType_fine)]
log_msg("Analysis groups: %s", paste(sort(unique(na.omit(md$AnalysisGroup))), collapse = ", "))

###############################################################################
# Step 2. Pseudobulk: sum raw counts per library x analysis group
###############################################################################
keep <- !is.na(md$AnalysisGroup)
pb_key <- paste(md$Run[keep], md$AnalysisGroup[keep], sep = "|")
counts <- LayerData(obj, assay = "RNA", layer = "counts")[, keep, drop = FALSE]

# Sparse matrix multiply is the cheap way to sum columns by group
design_pb <- sparse.model.matrix(~ 0 + factor(pb_key))
colnames(design_pb) <- levels(factor(pb_key))
pb <- as.matrix(counts %*% design_pb)
rm(counts); invisible(gc())

pb_meta <- data.frame(sample_id = colnames(pb), stringsAsFactors = FALSE) %>%
  separate(sample_id, into = c("Run", "AnalysisGroup"), sep = "\\|", remove = FALSE) %>%
  left_join(md %>% distinct(Run, Donor, Condition, Sex, LibraryType, Family,
                            Age = if ("Age" %in% colnames(md)) Age else NA_real_),
            by = "Run")
pb_meta$n_cells <- as.integer(table(pb_key)[pb_meta$sample_id])
write.csv(pb_meta, file.path(TF_DIR, "pseudobulk_summary.csv"), row.names = FALSE, na = "")
log_msg("Pseudobulk: %d samples x %d genes (%d groups)", ncol(pb), nrow(pb),
        n_distinct(pb_meta$AnalysisGroup))

###############################################################################
# Helper: fit one analysis group, return per-gene statistics
###############################################################################
# fast = TRUE skips duplicateCorrelation. Used only for the leave-one-out and
# permutation loops, where the same simplification applies to every refit and
# so does not bias the comparison. The reported fits always use fast = FALSE.
fit_group <- function(grp, block_var, perm_labels = NULL, drop_donor = NULL, fast = FALSE) {
  sel <- pb_meta$AnalysisGroup == grp & pb_meta$n_cells >= PB_MIN_CELLS
  if (!is.null(drop_donor)) sel <- sel & pb_meta$Donor != drop_donor
  m <- pb_meta[sel, , drop = FALSE]
  if (nrow(m) < 4) return(NULL)
  cond <- if (is.null(perm_labels)) m$Condition else perm_labels[m$Donor]
  if (n_distinct(cond) < 2) return(NULL)
  # both groups need at least 2 donors, else the contrast is meaningless
  if (min(table(tapply(cond, m$Donor, function(x) x[1]))) < 2) return(NULL)

  y <- DGEList(counts = pb[, m$sample_id, drop = FALSE])
  y <- y[filterByExpr(y, group = cond), , keep.lib.sizes = FALSE]
  if (nrow(y) < 500) return(NULL)
  y <- calcNormFactors(y)

  # Covariates are included only when they vary and are not aliased with condition
  cond_f <- factor(cond, levels = c("Healthy", "FPD"))
  terms <- list(Condition = cond_f)
  for (v in DE_COVARIATES) {
    if (!v %in% colnames(m)) next
    x <- m[[v]]
    if (all(is.na(x)) || n_distinct(na.omit(x)) < 2) next
    if (is.character(x) || is.factor(x)) {
      x <- factor(x)
      # skip a covariate that is perfectly confounded with condition
      if (nrow(unique(data.frame(x, cond_f))) <= n_distinct(x)) next
      terms[[v]] <- x
    } else {
      if (anyNA(x)) next
      terms[[v]] <- x
    }
  }
  design <- model.matrix(as.formula(paste("~", paste(names(terms), collapse = " + "))),
                         data = as.data.frame(terms))
  if (qr(design)$rank < ncol(design)) {          # drop aliased columns
    design <- design[, qr(design)$pivot[seq_len(qr(design)$rank)], drop = FALSE]
  }
  coef_name <- grep("^Condition", colnames(design), value = TRUE)[1]
  if (is.na(coef_name)) return(NULL)

  v <- voom(y, design)
  blk <- m[[block_var]]
  if (!fast && n_distinct(blk) < nrow(m) && n_distinct(blk) > 1) {
    dc <- suppressWarnings(duplicateCorrelation(v, design, block = blk))
    if (is.finite(dc$consensus.correlation) && abs(dc$consensus.correlation) < 0.99) {
      v  <- voom(y, design, block = blk, correlation = dc$consensus.correlation)
      dc <- suppressWarnings(duplicateCorrelation(v, design, block = blk))
      fit <- lmFit(v, design, block = blk, correlation = dc$consensus.correlation)
    } else fit <- lmFit(v, design)
  } else fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

  list(stats = data.frame(gene = rownames(tt), t = tt$t, logFC = tt$logFC,
                          P = tt$P.Value, adjP = tt$adj.P.Val, row.names = NULL),
       logcpm = cpm(y, log = TRUE), meta = m, coef = coef_name,
       n_FPD_donors = n_distinct(m$Donor[cond == "FPD"]),
       n_HD_donors  = n_distinct(m$Donor[cond == "Healthy"]),
       covariates   = paste(setdiff(names(terms), "Condition"), collapse = ";"))
}

# decoupleR ULM on a genes x contrasts matrix
run_activity <- function(mat, net, cols, minsize = DECOUPLER_MINSIZE) {
  run_ulm(mat = as.matrix(mat), network = net, .source = cols$source,
          .target = cols$target, .mor = cols$mor, minsize = minsize)
}

###############################################################################
# Step 3. Fit every group under both correlation structures
###############################################################################
grps <- sort(unique(na.omit(md$AnalysisGroup)))
results <- list(); perdonor <- list(); skipped <- list(); genestats <- list()

for (block_var in c("Donor", "Family")) {
  for (grp in grps) {
    fg <- fit_group(grp, block_var)
    if (is.null(fg)) {
      if (block_var == "Donor")
        skipped[[grp]] <- data.frame(AnalysisGroup = grp,
                                     reason = "too few samples/donors, or too few genes after filtering")
      next
    }
    log_msg("[block=%-6s] %-16s genes=%5d  FPD donors=%d  Healthy donors=%d  covars=[%s]",
            block_var, grp, nrow(fg$stats), fg$n_FPD_donors, fg$n_HD_donors, fg$covariates)

    # Per-gene statistics are kept: script 06 tests gene sets on these ranks
    genestats[[paste(block_var, grp)]] <- fg$stats %>%
      mutate(cohort = "BM_scRNA", cell_type = grp, block = block_var, .before = 1)

    tmat <- matrix(fg$stats$t, ncol = 1, dimnames = list(fg$stats$gene, grp))
    for (nm in c("CollecTRI", "PROGENy")) {
      net  <- if (nm == "CollecTRI") collectri else progeny
      cols <- if (nm == "CollecTRI") cc else pc
      act <- run_activity(tmat, net, cols)
      results[[paste(block_var, grp, nm)]] <- act %>%
        transmute(cohort = "BM_scRNA", stage = grp, cell_type = grp,
                  feature_type = if (nm == "CollecTRI") "TF_activity" else "pathway_activity",
                  feature = source, tf = source,
                  method = paste0("decoupleR_ULM_", nm), block = block_var,
                  effect = score, stat = score, p_value = p_value,
                  n_FPD = fg$n_FPD_donors, n_Healthy = fg$n_HD_donors,
                  covariates = fg$covariates)
    }

    # Per-donor scores, from each sample's logCPM (used for plots and for
    # correlating with cytokine levels later)
    if (block_var == "Donor") {
      act_s <- run_activity(fg$logcpm, collectri, cc)
      perdonor[[grp]] <- act_s %>%
        rename(sample_id = condition) %>%
        left_join(fg$meta %>% select(sample_id, Donor, Condition, Sex, Age, Family,
                                     LibraryType, n_cells), by = "sample_id") %>%
        transmute(cohort = "BM_scRNA", cell_type = grp, donor_id = Donor, family_id = Family,
                  condition = Condition, sex = Sex, age = Age, library_type = LibraryType,
                  n_cells, feature = source, score)
    }
  }
}

if (length(results) == 0) stop("No analysis group could be fitted. Check pseudobulk_summary.csv.")
res <- bind_rows(results)

###############################################################################
# Step 4. FDR in two families (pre-specified panel vs discovery)
###############################################################################
res$fdr_family <- ifelse(res$feature %in% TF_PANEL, "panel", "discovery")
res$fdr_family[res$feature_type == "pathway_activity"] <- "pathway"
res <- res %>% group_by(block, cell_type, feature_type, fdr_family) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>% ungroup()

# Panel TFs that could not be scored anywhere (too few targets in the network)
scored <- unique(res$feature[res$feature_type == "TF_activity"])
not_testable <- setdiff(TF_PANEL, scored)
if (length(not_testable) > 0) {
  log_msg("Panel TFs NOT testable (too few targets in CollecTRI, or not expressed): %s",
          paste(not_testable, collapse = ", "))
  res <- bind_rows(res, expand.grid(block = unique(res$block), cell_type = grps,
                                    feature = not_testable, stringsAsFactors = FALSE) %>%
    mutate(cohort = "BM_scRNA", stage = cell_type, feature_type = "TF_activity",
           tf = feature, method = "decoupleR_ULM_CollecTRI", effect = NA_real_,
           stat = NA_real_, p_value = NA_real_, p_adj = NA_real_,
           fdr_family = "panel", notes = "not testable: too few targets"))
}
res$network_version <- NETWORK_VERSION_TAG

###############################################################################
# Step 5. Leave-one-donor-out sign stability (panel TFs only; it is slow)
###############################################################################
log_msg("Leave-one-donor-out stability for %d panel TFs", length(intersect(TF_PANEL, scored)))
loo <- list()
for (grp in grps) {
  base <- res %>% filter(block == "Donor", cell_type == grp, feature_type == "TF_activity",
                         fdr_family == "panel", !is.na(effect))
  if (nrow(base) == 0) next
  donors <- unique(pb_meta$Donor[pb_meta$AnalysisGroup == grp])
  signs <- list()
  for (d in donors) {
    fg <- tryCatch(fit_group(grp, "Donor", drop_donor = d, fast = TRUE), error = function(e) NULL)
    if (is.null(fg)) next
    tmat <- matrix(fg$stats$t, ncol = 1, dimnames = list(fg$stats$gene, grp))
    a <- run_activity(tmat, collectri, cc)
    signs[[d]] <- setNames(sign(a$score), a$source)
  }
  if (length(signs) == 0) next
  stable <- sapply(base$feature, function(f) {
    s <- sapply(signs, function(v) if (f %in% names(v)) v[[f]] else NA_real_)
    all(na.omit(s) == sign(base$effect[base$feature == f][1])) && sum(!is.na(s)) > 0
  })
  loo[[grp]] <- data.frame(cell_type = grp, feature = base$feature, loo_sign_stable = stable)
}
loo <- bind_rows(loo)
res <- res %>% left_join(loo, by = c("cell_type", "feature"))

###############################################################################
# Step 6. Permutation null (condition labels shuffled among donors, families kept whole)
###############################################################################
log_msg("Permutations: %d, shuffling condition among families", N_PERM)
donor_tab <- pb_meta %>% distinct(Donor, Family, Condition)
fam_tab   <- donor_tab %>% distinct(Family, Condition)
perm_rows <- list()
set.seed(SEED)
for (i in seq_len(N_PERM)) {
  shuffled <- fam_tab
  shuffled$Condition <- sample(shuffled$Condition)         # whole families move together
  lab <- setNames(shuffled$Condition[match(donor_tab$Family, shuffled$Family)], donor_tab$Donor)
  for (grp in grps) {
    fg <- tryCatch(fit_group(grp, "Donor", perm_labels = lab, fast = TRUE), error = function(e) NULL)
    if (is.null(fg)) next
    tmat <- matrix(fg$stats$t, ncol = 1, dimnames = list(fg$stats$gene, grp))
    a <- run_activity(tmat, collectri, cc)
    perm_rows[[length(perm_rows) + 1]] <- data.frame(perm = i, cell_type = grp,
                                                     feature = a$source, effect = a$score)
  }
  if (i %% 25 == 0) log_msg("  permutation %d/%d", i, N_PERM)
}
perm <- bind_rows(perm_rows)

###############################################################################
# Step 7. Write the standard tables
###############################################################################
stem <- file.path(TF_DIR, "BM_scRNA_progenitors_decoupleR")
write.table(res,      paste0(stem, "_results.tsv"),  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(perdonor), paste0(stem, "_perdonor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(perm) > 0)
  write.table(perm,   paste0(stem, "_perm.tsv"),     sep = "\t", quote = FALSE, row.names = FALSE)
if (length(skipped) > 0)
  write.csv(bind_rows(skipped), file.path(TF_DIR, "groups_skipped.csv"), row.names = FALSE)
write.table(bind_rows(genestats), file.path(TF_DIR, "per_gene_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Wrote %d result rows, %d per-donor rows, %d permutation rows",
        nrow(res), nrow(bind_rows(perdonor)), nrow(perm))

###############################################################################
# Step 8. Plots
###############################################################################
panel_res <- res %>% filter(feature_type == "TF_activity", fdr_family == "panel", !is.na(effect))
if (nrow(panel_res) > 0) {
  p_heat <- ggplot(panel_res, aes(cell_type, factor(feature, levels = rev(TF_PANEL)), fill = effect)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(p_adj < 0.05, "*", "")), vjust = 0.75, size = 4) +
    facet_wrap(~ block) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                         name = "ULM score\n(+ = higher in FPD)") +
    labs(x = NULL, y = NULL, title = "Pre-specified TF panel: activity in FPD vs Healthy",
         subtitle = "* = adjusted p < 0.05 within the panel family; the two panels are the two correlation structures") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_plot(p_heat, PLOT_DIR, "panel_TF_activity_heatmap", width = 10,
            height = max(5, 0.22 * length(TF_PANEL) + 2))
}

path_res <- res %>% filter(feature_type == "pathway_activity", !is.na(effect))
if (nrow(path_res) > 0) {
  p_path <- ggplot(path_res, aes(cell_type, feature, fill = effect)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(p_adj < 0.05, "*", "")), vjust = 0.75, size = 4) +
    facet_wrap(~ block) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                         name = "ULM score") +
    labs(x = NULL, y = NULL, title = "PROGENy pathway activity, FPD vs Healthy") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_plot(p_path, PLOT_DIR, "PROGENy_pathway_activity_heatmap", width = 10, height = 6)
}

# Per-donor dots for the panel TFs that move most, so single-donor effects show
pd <- bind_rows(perdonor)
if (nrow(pd) > 0 && nrow(panel_res) > 0) {
  top_tf <- panel_res %>% filter(block == "Donor") %>% slice_min(p_value, n = 12) %>% pull(feature)
  pdd <- pd %>% filter(feature %in% top_tf)
  if (nrow(pdd) > 0) {
    p_dot <- ggplot(pdd, aes(condition, score, colour = condition)) +
      geom_point(position = position_jitter(width = 0.15, seed = SEED), size = 1.8) +
      facet_grid(feature ~ cell_type, scales = "free_y") +
      scale_colour_manual(values = CONDITION_COLORS) +
      labs(x = NULL, y = "Per-sample TF activity (ULM on logCPM)",
           title = "Per-donor TF activity for the strongest panel hits",
           subtitle = "One point per library. Check that an effect is not carried by a single donor.") +
      theme_bw(base_size = 8) + theme(legend.position = "none")
    save_plot(p_dot, PLOT_DIR, "panel_TF_per_donor", width = 12, height = 14, pdf = FALSE)
  }
}

log_mem("end of 05")
log_disk()
write_session_info("05_tf_activity")
log_msg("=== 05_tf_activity.R finished ===")
