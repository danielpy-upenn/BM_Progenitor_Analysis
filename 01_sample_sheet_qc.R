###############################################################################
# 01_sample_sheet_qc.R
#
# What this script does
#   1. Reads SraRunTable.csv and builds a sample sheet:
#        Run -> Donor (submitted_subject_id) -> Condition (subject_is_affected)
#        -> Sex, plus Sample.Name.
#   2. Reads metadata/run_library_map.csv (filled in by you) to link each run
#      to its library in the original preprocessing.R, which sets its QC
#      thresholds and library type. On the FIRST run this file does not exist
#      yet: the script writes a template and stops.
#   3. Loads each run's Cell Ranger filtered matrix into Seurat, computes
#      percent.mt, draws the QC scatter plots with the thresholds drawn,
#      filters cells, and records cells before/after.
#   4. Merges all runs into one Seurat v5 object with RNA layers split by run,
#      and saves a checkpoint for script 02.
#
# Inputs : SRR_Folders/SRR*/outs/filtered_feature_bc_matrix/, SraRunTable.csv,
#          metadata/qc_thresholds.csv, metadata/run_library_map.csv
# Outputs: Output/01_qc/ (sample_sheet.csv, qc_summary.csv, per-run QC plots)
#          Output/checkpoints/01_qc_filtered_merged.rds
#
# Run    : Rscript 01_sample_sheet_qc.R
###############################################################################

# ---- Load config (00_config.R sits in the same folder as this script) --------
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
  library(patchwork)
})

ensure_dirs(QC_DIR, file.path(QC_DIR, "per_run_plots"), LOG_DIR, CKPT_DIR, META_DIR)
log_msg("=== 01_sample_sheet_qc.R started ===")
log_msg("Seurat version: %s", as.character(packageVersion("Seurat")))
gc(reset = TRUE)

###############################################################################
# Step 1. Sample sheet from SraRunTable.csv
###############################################################################
if (!file.exists(SRA_TABLE)) stop("SraRunTable.csv not found at ", SRA_TABLE)

# read.csv turns the header "Sample Name" into "Sample.Name"
sra <- read.csv(SRA_TABLE, stringsAsFactors = FALSE, check.names = TRUE)

required_cols <- c("Run", "submitted_subject_id", "subject_is_affected", "sex",
                   "Sample.Name")
missing_cols <- setdiff(required_cols, colnames(sra))
if (length(missing_cols) > 0) {
  stop("SraRunTable.csv is missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nColumns found: ", paste(colnames(sra), collapse = ", "))
}

sra <- sra %>% mutate(across(all_of(required_cols), ~ trimws(as.character(.x))))

# Keep only the 30 runs of this study and check none are missing
missing_runs <- setdiff(EXPECTED_RUNS, sra$Run)
if (length(missing_runs) > 0) {
  stop("These expected runs are not in SraRunTable.csv: ", paste(missing_runs, collapse = ", "))
}
sra <- sra %>% filter(Run %in% EXPECTED_RUNS)
if (any(duplicated(sra$Run))) stop("Duplicated Run IDs in SraRunTable.csv")

# Fix known typos in the subject ID (FD19006 -> HD19006), with a safety check
for (wrong in names(SUBJECT_ID_FIXES)) {
  idx <- sra$submitted_subject_id == wrong
  if (any(idx)) {
    right <- SUBJECT_ID_FIXES[[wrong]]
    if (grepl("^HD", right) && any(tolower(sra$subject_is_affected[idx]) != "no")) {
      stop(sprintf("Refusing to relabel %s as %s: subject_is_affected is not 'No'.", wrong, right))
    }
    sra$submitted_subject_id[idx] <- right
    log_msg("Relabeled submitted_subject_id %s -> %s for %d run(s)", wrong, right, sum(idx))
  }
}

# Recode condition
affected <- tolower(sra$subject_is_affected)
if (!all(affected %in% c("yes", "no"))) {
  stop("Unexpected values in subject_is_affected: ",
       paste(unique(sra$subject_is_affected[!affected %in% c("yes", "no")]), collapse = ", "))
}

sample_sheet <- sra %>%
  transmute(Run,
            Donor       = submitted_subject_id,
            Condition   = ifelse(tolower(subject_is_affected) == "yes", "FPD", "Healthy"),
            Sex         = tolower(sex),
            Sample.Name = Sample.Name) %>%
  arrange(Run)

# Each donor must have exactly one condition and one sex
donor_check <- sample_sheet %>%
  group_by(Donor) %>%
  summarise(n_cond = n_distinct(Condition), n_sex = n_distinct(Sex),
            Condition = first(Condition), n_runs = n(), .groups = "drop")
bad <- donor_check %>% filter(n_cond > 1 | n_sex > 1)
if (nrow(bad) > 0) {
  stop("Donor(s) with conflicting condition or sex across runs: ", paste(bad$Donor, collapse = ", "))
}
n_fpd <- sum(donor_check$Condition == "FPD")
n_hd  <- sum(donor_check$Condition == "Healthy")
log_msg("Donors: %d FPD, %d Healthy (expected %d and %d)", n_fpd, n_hd,
        EXPECTED_N_FPD, EXPECTED_N_HEALTHY)
if (n_fpd != EXPECTED_N_FPD || n_hd != EXPECTED_N_HEALTHY) {
  stop("Donor counts do not match the study design. Check submitted_subject_id / subject_is_affected.")
}
print(as.data.frame(donor_check %>% select(Donor, Condition, n_runs)))

# Every run folder must exist
sample_sheet$matrix_dir <- file.path(SRR_DIR, sample_sheet$Run, "outs", "filtered_feature_bc_matrix")
missing_dirs <- sample_sheet$Run[!dir.exists(sample_sheet$matrix_dir)]
if (length(missing_dirs) > 0) {
  stop("filtered_feature_bc_matrix folder not found for: ", paste(missing_dirs, collapse = ", "))
}

###############################################################################
# Step 2. Link runs to libraries (QC thresholds + library type)
###############################################################################
qc_thr <- read.csv(QC_THRESHOLDS_CSV, stringsAsFactors = FALSE)
qc_thr$mt_min <- suppressWarnings(as.numeric(qc_thr$mt_min))

if (!file.exists(RUN_LIBRARY_MAP_CSV)) {
  template <- sample_sheet %>%
    transmute(Run, Donor, Condition, Sample.Name, Library = "", LibraryType = "")
  write.csv(template, RUN_LIBRARY_TEMPLATE, row.names = FALSE, na = "")
  log_msg("Wrote %s", RUN_LIBRARY_TEMPLATE)
  stop(paste0(
    "\n\nrun_library_map.csv does not exist yet. A template was written to:\n  ",
    RUN_LIBRARY_TEMPLATE, "\n",
    "For each Run, fill in:\n",
    "  Library     = the matching label from the 'Library' column of qc_thresholds.csv\n",
    "                (leave blank to use the default thresholds)\n",
    "  LibraryType = one of: ", paste(LIBRARY_TYPES, collapse = ", "), "\n",
    "                (required when Library is blank; optional otherwise)\n",
    "Save it as ", RUN_LIBRARY_MAP_CSV, " and rerun this script.\n"))
}

run_map <- read.csv(RUN_LIBRARY_MAP_CSV, stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (!all(c("Run", "Library", "LibraryType") %in% colnames(run_map))) {
  stop("run_library_map.csv must have columns Run, Library, LibraryType")
}
run_map <- run_map %>%
  mutate(Run = trimws(Run), Library = trimws(Library), LibraryType = trimws(LibraryType)) %>%
  select(Run, Library, LibraryType)

# --- Validation of the map ---
problems <- character(0)
if (!setequal(run_map$Run, sample_sheet$Run) || any(duplicated(run_map$Run))) {
  problems <- c(problems, sprintf(
    "Runs in map and run table differ or are duplicated. Missing: [%s]  Extra: [%s]",
    paste(setdiff(sample_sheet$Run, run_map$Run), collapse = ", "),
    paste(setdiff(run_map$Run, sample_sheet$Run), collapse = ", ")))
}
lib_given <- run_map$Library[!is.na(run_map$Library)]
unknown_libs <- setdiff(lib_given, qc_thr$Library)
if (length(unknown_libs) > 0) {
  problems <- c(problems, paste("Library label(s) not found in qc_thresholds.csv (typo?):",
                                paste(unknown_libs, collapse = ", ")))
}
dup_libs <- unique(lib_given[duplicated(lib_given)])
if (length(dup_libs) > 0) {
  problems <- c(problems, paste("Library label(s) used for more than one run:",
                                paste(dup_libs, collapse = ", ")))
}
bad_types <- setdiff(na.omit(run_map$LibraryType), LIBRARY_TYPES)
if (length(bad_types) > 0) {
  problems <- c(problems, paste("Unknown LibraryType value(s):", paste(bad_types, collapse = ", "),
                                "| allowed:", paste(LIBRARY_TYPES, collapse = ", ")))
}
no_type <- run_map$Run[is.na(run_map$Library) & is.na(run_map$LibraryType)]
if (length(no_type) > 0) {
  problems <- c(problems, paste("Run(s) with blank Library must have a LibraryType:",
                                paste(no_type, collapse = ", ")))
}
# The library's donor must match the run's donor (catches mis-assignments)
chk <- sample_sheet %>% select(Run, Donor) %>%
  inner_join(run_map, by = "Run") %>%
  filter(!is.na(Library)) %>%
  left_join(qc_thr %>% select(Library, LibDonor = Donor), by = "Library") %>%
  filter(!is.na(LibDonor) & LibDonor != Donor)
if (nrow(chk) > 0) {
  problems <- c(problems, paste0("Library assigned to the wrong donor: ",
                                 paste(sprintf("%s (donor %s) <- %s (donor %s)", chk$Run, chk$Donor,
                                               chk$Library, chk$LibDonor), collapse = "; ")))
}
if (length(problems) > 0) {
  stop("run_library_map.csv has problems:\n  - ", paste(problems, collapse = "\n  - "))
}

# --- Build final per-run QC table ---
sample_sheet <- sample_sheet %>%
  left_join(run_map, by = "Run") %>%
  left_join(qc_thr %>% select(Library, LibTypeTable = LibraryType, nFeature_min,
                              nFeature_max, mt_max, mt_min), by = "Library")

# Library type: your map wins if filled in, otherwise the thresholds table
type_conflict <- !is.na(sample_sheet$LibraryType) & !is.na(sample_sheet$LibTypeTable) &
  sample_sheet$LibraryType != sample_sheet$LibTypeTable
if (any(type_conflict)) {
  log_msg("NOTE: LibraryType in run_library_map.csv overrides qc_thresholds.csv for: %s",
          paste(sample_sheet$Run[type_conflict], collapse = ", "))
}
sample_sheet$LibraryType <- ifelse(is.na(sample_sheet$LibraryType),
                                   sample_sheet$LibTypeTable, sample_sheet$LibraryType)
sample_sheet$LibTypeTable <- NULL

# Defaults for runs with no library match
use_default <- is.na(sample_sheet$Library)
sample_sheet$QC_source <- ifelse(use_default, "default", "preprocessing.R")
sample_sheet$nFeature_min[use_default] <- DEFAULT_QC$nFeature_min
sample_sheet$nFeature_max[use_default] <- DEFAULT_QC$nFeature_max
sample_sheet$mt_max[use_default]       <- DEFAULT_QC$mt_max
sample_sheet$mt_min[use_default]       <- DEFAULT_QC$mt_min
sample_sheet$Library[use_default]      <- paste0(sample_sheet$Donor[use_default], "_",
                                                 sample_sheet$Run[use_default])
if (any(use_default)) {
  log_msg("Runs using DEFAULT thresholds (nFeature %d-%d, percent.mt < %d): %s",
          DEFAULT_QC$nFeature_min, DEFAULT_QC$nFeature_max, DEFAULT_QC$mt_max,
          paste(sample_sheet$Run[use_default], collapse = ", "))
}

write.csv(sample_sheet %>% select(-matrix_dir), file.path(QC_DIR, "sample_sheet.csv"),
          row.names = FALSE, na = "")
log_msg("Sample sheet written (%d runs)", nrow(sample_sheet))

###############################################################################
# Step 3. Load each run, QC plots, filter
###############################################################################
qc_scatter <- function(obj, thr, title) {
  md <- obj@meta.data
  mt_lines <- c(thr$mt_max, if (!is.na(thr$mt_min)) thr$mt_min)
  p1 <- ggplot(md, aes(nCount_RNA, percent.mt)) +
    geom_point(size = 0.3, alpha = 0.4, colour = "grey25") +
    geom_hline(yintercept = mt_lines, linetype = "dashed", colour = "#C0392B") +
    labs(x = "nCount_RNA", y = "percent.mt") + theme_classic()
  p2 <- ggplot(md, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(size = 0.3, alpha = 0.4, colour = "grey25") +
    geom_hline(yintercept = c(thr$nFeature_min, thr$nFeature_max),
               linetype = "dashed", colour = "#C0392B") +
    labs(x = "nCount_RNA", y = "nFeature_RNA") + theme_classic()
  (p1 | p2) + plot_annotation(title = title,
                              subtitle = "Before filtering; dashed lines = thresholds")
}

obj_list   <- vector("list", nrow(sample_sheet))
qc_summary <- vector("list", nrow(sample_sheet))

for (i in seq_len(nrow(sample_sheet))) {
  s <- sample_sheet[i, ]
  log_msg("[%02d/%02d] %s | %s | %s | %s", i, nrow(sample_sheet), s$Run, s$Donor,
          s$LibraryType, s$QC_source)

  counts <- Read10X(data.dir = s$matrix_dir)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]   # multi-feature output
  colnames(counts) <- paste(s$Run, colnames(counts), sep = "_") # unique barcodes

  obj <- CreateSeuratObject(counts = counts, project = s$Run,
                            min.cells = MIN_CELLS_PER_GENE,
                            min.features = MIN_FEATURES_PER_CELL)
  rm(counts)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
  n_mt_genes <- sum(grepl(MT_PATTERN, rownames(obj)))
  if (n_mt_genes == 0) warning(s$Run, ": no genes match ", MT_PATTERN)

  # Metadata
  for (col in c("Run", "Donor", "Condition", "Sex", "Library", "LibraryType",
                "Sample.Name", "QC_source")) {
    obj[[col]] <- s[[col]]
  }

  # QC plot before filtering
  title <- sprintf("%s | %s | %s | %s", s$Run, s$Donor, s$Library, s$LibraryType)
  save_plot(qc_scatter(obj, s, title), file.path(QC_DIR, "per_run_plots"),
            paste0(s$Run, "_QC_scatter"), width = 10, height = 4.5, pdf = FALSE)

  n_before <- ncol(obj)
  md <- obj@meta.data
  keep <- md$nFeature_RNA > s$nFeature_min & md$nFeature_RNA < s$nFeature_max &
    md$percent.mt < s$mt_max
  if (!is.na(s$mt_min)) keep <- keep & md$percent.mt > s$mt_min
  if (sum(keep) == 0) stop(s$Run, ": zero cells left after filtering; check thresholds")
  obj <- subset(obj, cells = colnames(obj)[keep])

  qc_summary[[i]] <- data.frame(
    Run = s$Run, Donor = s$Donor, Condition = s$Condition, Sex = s$Sex,
    Library = s$Library, LibraryType = s$LibraryType, QC_source = s$QC_source,
    nFeature_min = s$nFeature_min, nFeature_max = s$nFeature_max,
    mt_max = s$mt_max, mt_min = s$mt_min,
    cells_before = n_before, cells_after = ncol(obj),
    pct_kept = round(100 * ncol(obj) / n_before, 1),
    median_nFeature_after = median(obj$nFeature_RNA),
    median_nCount_after = median(obj$nCount_RNA),
    median_pct_mt_after = round(median(obj$percent.mt), 2))
  log_msg("    cells: %d -> %d (%.1f%% kept)", n_before, ncol(obj), 100 * ncol(obj) / n_before)
  obj_list[[i]] <- obj
}

qc_summary <- bind_rows(qc_summary)
write.csv(qc_summary, file.path(QC_DIR, "qc_summary.csv"), row.names = FALSE, na = "")
log_msg("Total cells after QC: %d", sum(qc_summary$cells_after))

# ---- Summary plots -----------------------------------------------------------
p_cells_run <- ggplot(qc_summary, aes(x = reorder(Run, cells_after), y = cells_after,
                                      fill = Condition)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = paste(Donor, LibraryType)), hjust = -0.05, size = 2.6) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = NULL, y = "Cells after QC", title = "Cells per run after QC") +
  theme_classic()
save_plot(p_cells_run, QC_DIR, "cells_per_run_after_QC", width = 8, height = 8)

donor_tot <- qc_summary %>% group_by(Donor, Condition, Sex) %>%
  summarise(cells = sum(cells_after), .groups = "drop")
p_cells_donor <- ggplot(donor_tot, aes(x = reorder(Donor, cells), y = cells, fill = Condition)) +
  geom_col(width = 0.8) + coord_flip() +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Cells after QC (all libraries)", title = "Cells per donor after QC") +
  theme_classic()
save_plot(p_cells_donor, QC_DIR, "cells_per_donor_after_QC", width = 7, height = 5)

###############################################################################
# Step 4. Merge and split RNA layers by run
###############################################################################
log_msg("Merging %d runs", length(obj_list))
merged <- merge(x = obj_list[[1]], y = obj_list[-1], merge.data = FALSE)
rm(obj_list); invisible(gc())

# Make the layer layout deterministic: one counts layer per run, named by SRR
merged <- JoinLayers(merged)
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$Run)
log_msg("Merged object: %d genes x %d cells; layers: %s", nrow(merged), ncol(merged),
        paste(head(Layers(merged), 3), collapse = ", "))

# Violin summary across runs (after filtering)
p_vln <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "Run", pt.size = 0, ncol = 1, layer = "counts") &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
save_plot(p_vln, QC_DIR, "QC_violins_after_filtering_by_run", width = 12, height = 11)

saveRDS(merged, CKPT_01)
log_msg("Saved checkpoint %s", CKPT_01)
log_mem("end of 01")
log_disk()
write_session_info("01_sample_sheet_qc")
log_msg("=== 01_sample_sheet_qc.R finished ===")
