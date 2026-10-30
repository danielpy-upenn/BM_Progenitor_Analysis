###############################################################################
# 06_heptad_enrichment.R
#
# Are genes that RUNX1 occupies in human CD34+ HSPCs dysregulated in FPD?
#
# Binding comes from Beck et al., Blood 2013 (10.1182/blood-2013-03-490425):
# 28,793 regions in primary human CD34+ HSPCs with binding calls for seven TFs
# (RUNX1, FLI1, ERG, GATA2, TAL1/SCL, LYL1, LMO2).
#
# Three gene sets are tested, in increasing specificity:
#   RUNX1_bound       9,727 genes (75% of the table's genes)
#   RUNX1_ETS_bound   6,877 genes (RUNX1 with FLI1 and/or ERG: composite sites)
#   Heptad_all7         358 genes (all seven factors)
# The GRADIENT across the three is the result. RUNX1_bound alone covers three
# quarters of the background, so on its own it can barely discriminate.
#
# Method: limma::cameraPR, a rank-based competitive test on the per-gene
# t-statistic. No fold-change cutoff, and it corrects for inter-gene
# correlation, which matters a great deal for a set of 9,727 genes.
# Background = the genes in the Beck table that are testable in this dataset,
# never the whole genome.
#
# Also writes the binding regions lifted to hg38, for the ATAC analysis.
#
# Inputs : Output/05_tf_activity/ (per-gene statistics), reference/beck_*
# Outputs: Output/06_heptad/
# Run    : Rscript 06_heptad_enrichment.R
###############################################################################

get_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}
source(file.path(get_script_dir(), "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(limma); library(ggplot2)
})

HEP_DIR <- file.path(OUT_DIR, "06_heptad")
ensure_dirs(HEP_DIR, LOG_DIR)
log_msg("=== 06_heptad_enrichment.R started ===")

###############################################################################
# Step 0. Inputs
###############################################################################
if (!file.exists(HEPTAD_GENESETS)) stop("Missing ", HEPTAD_GENESETS)
hept <- read_tsv(HEPTAD_GENESETS, show_col_types = FALSE)
missing_sets <- setdiff(HEPTAD_SETS, colnames(hept))
if (length(missing_sets) > 0) stop("Gene-set columns missing: ", paste(missing_sets, collapse = ", "))
log_msg("Heptad table: %d genes", nrow(hept))
for (s in HEPTAD_SETS) log_msg("  %-18s %5d genes", s, sum(hept[[s]] == 1))

# Per-gene statistics from script 05. That script does not currently persist
# them, so they are recomputed here from the same pseudobulk fit if absent.
GENESTAT_FILE <- file.path(OUT_DIR, "05_tf_activity", "per_gene_statistics.tsv")
if (!file.exists(GENESTAT_FILE)) {
  stop("Per-gene statistics not found: ", GENESTAT_FILE,
       "\nRun 05_tf_activity.R first (it writes this file).")
}
gs <- read_tsv(GENESTAT_FILE, show_col_types = FALSE)
req <- c("cell_type", "block", "gene", "t")
if (!all(req %in% colnames(gs))) stop("per_gene_statistics.tsv needs: ", paste(req, collapse = ", "))
log_msg("Per-gene statistics: %d rows, %d cell types",
        nrow(gs), n_distinct(gs$cell_type))

###############################################################################
# Step 1. cameraPR per cell type x gene set
###############################################################################
run_one <- function(df, set_name) {
  # Background: genes present in BOTH the Beck table and this cell type's fit
  bg <- intersect(df$gene, hept$gene)
  if (length(bg) < 200) return(NULL)
  stat <- setNames(df$t[match(bg, df$gene)], bg)
  in_set <- hept$gene[hept[[set_name]] == 1]
  idx <- which(bg %in% in_set)
  if (length(idx) < 10 || length(idx) > length(bg) - 10) return(NULL)
  out <- cameraPR(stat, list(set = idx), use.ranks = FALSE)
  data.frame(gene_set = set_name,
             n_background = length(bg), n_in_set = length(idx),
             frac_of_background = round(length(idx) / length(bg), 3),
             Direction = out$Direction[1], p_value = out$PValue[1],
             mean_t_in_set = mean(stat[idx]), mean_t_out = mean(stat[-idx]),
             row.names = NULL)
}

res <- list()
for (blk in unique(gs$block)) {
  for (ct in unique(gs$cell_type)) {
    df <- gs %>% filter(block == blk, cell_type == ct, is.finite(t))
    if (nrow(df) < 500) next
    for (s in HEPTAD_SETS) {
      r <- run_one(df, s)
      if (!is.null(r)) res[[paste(blk, ct, s)]] <-
        cbind(data.frame(cohort = "BM_scRNA", cell_type = ct, block = blk), r)
    }
  }
}
if (length(res) == 0) stop("No cell type produced a testable gene set.")
res <- bind_rows(res) %>%
  group_by(block, cell_type) %>% mutate(p_adj = p.adjust(p_value, "BH")) %>% ungroup() %>%
  mutate(method = "limma_cameraPR",
         network_version = "Beck2013_heptad_CD34_hg19",
         signed_effect = mean_t_in_set - mean_t_out) %>%
  arrange(block, cell_type, match(gene_set, HEPTAD_SETS))

write.table(res, file.path(HEP_DIR, "BM_scRNA_heptad_cameraPR_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(res %>% filter(block == "Donor") %>%
        select(cell_type, gene_set, n_in_set, Direction, p_value, p_adj, signed_effect)))

###############################################################################
# Step 2. Plot the gradient. A real RUNX1 effect should strengthen from the
# broad set to the all-seven set; a flat line across the three is a warning
# that the signal is not specific to RUNX1 occupancy.
###############################################################################
res$gene_set <- factor(res$gene_set, levels = HEPTAD_SETS)
p <- ggplot(res, aes(gene_set, signed_effect, group = cell_type, colour = cell_type)) +
  geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(aes(size = ifelse(p_adj < 0.05, "adj p < 0.05", "n.s."))) +
  scale_size_manual(values = c("adj p < 0.05" = 3.2, "n.s." = 1.6), name = NULL) +
  facet_wrap(~ block) +
  labs(x = NULL, y = "mean t in set  -  mean t outside set",
       title = "Heptad binding sets: are RUNX1-occupied genes dysregulated in FPD?",
       subtitle = "Specificity increases left to right. A rising line supports a RUNX1-occupancy effect.") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p, HEP_DIR, "heptad_set_gradient", width = 10, height = 5)

###############################################################################
# Step 3. Lift the binding regions to hg38 for the ATAC analysis
#
# The Beck coordinates are hg19 (verified: every chromosome's maximum end
# coordinate fits hg19 and exceeds hg38). ATAC peaks called on hg38 cannot be
# overlapped with them directly.
###############################################################################
if (!file.exists(HEPTAD_REGIONS)) {
  log_msg("Region BED not found (%s); skipping liftover", HEPTAD_REGIONS)
} else if (!file.exists(LIFTOVER_CHAIN)) {
  log_msg("Liftover chain not found at %s; skipping.", LIFTOVER_CHAIN)
  log_msg("  Download hg19ToHg38.over.chain.gz from UCSC, gunzip it, and place it there.")
} else if (!requireNamespace("rtracklayer", quietly = TRUE) ||
           !requireNamespace("GenomicRanges", quietly = TRUE)) {
  log_msg("rtracklayer/GenomicRanges not installed; skipping liftover")
} else {
  suppressPackageStartupMessages({ library(rtracklayer); library(GenomicRanges) })
  bed <- read_tsv(HEPTAD_REGIONS, show_col_types = FALSE)
  gr <- GRanges(bed$chr, IRanges(bed$start + 1L, bed$end))   # BED is 0-based, GRanges 1-based
  mcols(gr) <- bed[, setdiff(colnames(bed), c("chr", "start", "end"))]
  ch <- import.chain(LIFTOVER_CHAIN)
  lifted <- unlist(liftOver(gr, ch))
  # A region that maps to several places is ambiguous; keep only 1-to-1 mappings
  tab <- table(lifted$region_id)
  lifted <- lifted[lifted$region_id %in% names(tab)[tab == 1]]
  out <- data.frame(chr = as.character(seqnames(lifted)),
                    start = start(lifted) - 1L, end = end(lifted),
                    as.data.frame(mcols(lifted)))
  write.table(out, file.path(HEP_DIR, "beck_heptad_regions_hg38.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  log_msg("Liftover hg19 -> hg38: %d of %d regions mapped uniquely (%.1f%%)",
          nrow(out), nrow(bed), 100 * nrow(out) / nrow(bed))
  log_msg("  Give Output/06_heptad/beck_heptad_regions_hg38.bed to the ATAC analysis.")
}

write_session_info("06_heptad_enrichment")
log_msg("=== 06_heptad_enrichment.R finished ===")
