# Reference per-sample signature scores from singscore::simpleScore on a
# planted fixture. Atman's `score signatures --method singscore` must
# reproduce the TotalScore values within tight numerical tolerance.
#
# Run: Rscript singscore_reference.R
# (Set ATMAN_FIXTURE_DIR to override the output directory.)

suppressPackageStartupMessages({
  library(singscore)
})

script_arg <- Sys.getenv("ATMAN_FIXTURE_DIR")
fixture_dir <- if (nzchar(script_arg)) script_arg else {
  args <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(self) == 1 && file.exists(self)) dirname(normalizePath(self))
  else "crates/atman/tests/fixtures"
}

set.seed(20260428L)

# 4 samples × 30 genes. Three signatures planted:
#   * "responder_up" — genes that should rank high (positive score)
#   * "responder_down" — genes that should rank low (negative score)
#   * "scattered" — random subset (score near 0)
n_samples <- 4L
n_genes <- 30L
sample_ids <- sprintf("S%02d", seq_len(n_samples))
gene_ids <- sprintf("G%02d", seq_len(n_genes))

# Construct a sample × gene matrix where the per-sample ranking varies.
# The first 5 genes are "down" (low across all samples),
# the last 5 genes are "up" (high across all samples), middle is noise.
mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)
rownames(mat) <- gene_ids
colnames(mat) <- sample_ids
mat[1:5, ]            <- mat[1:5, ]  - 4   # down genes
mat[(n_genes - 4):n_genes, ] <- mat[(n_genes - 4):n_genes, ] + 4   # up genes

up_sig    <- sprintf("G%02d", (n_genes - 4):n_genes)
down_sig  <- sprintf("G%02d", 1:5)
scattered <- sprintf("G%02d", seq(8L, n_genes - 6L, by = 4L))  # ~6 mid genes

# Emit canonical Atman fixture files.
samples_tsv <- data.frame(
  sample_id = sample_ids,
  subject_id = sample_ids,
  condition = "Case",
  is_control = 0L,
  sample_type = "csf",
  ingest_order = seq_len(n_samples),
  stringsAsFactors = FALSE
)
write.table(samples_tsv, file = file.path(fixture_dir, "singscore_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

proteins_tsv <- data.frame(
  platform = "spectronaut_report",
  assay_id = sprintf("P%05d", seq_len(n_genes)),
  uniprot = sprintf("P%05d", seq_len(n_genes)),
  gene_symbol = gene_ids,
  panel = "spectronaut",
  panel_lot = "",
  stringsAsFactors = FALSE
)
write.table(proteins_tsv, file = file.path(fixture_dir, "singscore_proteins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Long-format measurements.tsv. One row per (sample, gene).
rows <- vector("list", n_samples * n_genes)
order_counter <- 0L
for (s_idx in seq_len(n_samples)) {
  for (g_idx in seq_len(n_genes)) {
    order_counter <- order_counter + 1L
    val <- sprintf("%.6f", mat[g_idx, s_idx])
    rows[[order_counter]] <- data.frame(
      platform = "spectronaut_report",
      sample_id = sample_ids[s_idx],
      assay_id = sprintf("P%05d", g_idx),
      gene_symbol = gene_ids[g_idx],
      panel = "spectronaut",
      npx_source_str = val,
      abundance = val,
      abundance_raw = val,
      abundance_unit = "log2_intensity",
      qc_sample = "PASS",
      qc_assay = "PASS",
      detection_limit = "",
      below_lod = 0L,
      dropped_by_qc = 0L,
      plate_id = "",
      panel_lot = "",
      ingest_order = order_counter,
      stringsAsFactors = FALSE
    )
  }
}
measurements <- do.call(rbind, rows)
write.table(measurements,
            file = file.path(fixture_dir, "singscore_measurements.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Long gene_sets TSV (set_name, gene_symbol).
sets_long <- rbind(
  data.frame(set_name = "responder_up",   gene_symbol = up_sig,    stringsAsFactors = FALSE),
  data.frame(set_name = "responder_down", gene_symbol = down_sig,  stringsAsFactors = FALSE),
  data.frame(set_name = "scattered",      gene_symbol = scattered, stringsAsFactors = FALSE)
)
write.table(sets_long,
            file = file.path(fixture_dir, "singscore_gene_sets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Compute the singscore reference. simpleScore with knownDirection=TRUE
# and centerScore=TRUE (defaults): TotalScore in [-0.5, 0.5].
ranked <- rankGenes(mat)
score_one <- function(genes) {
  res <- simpleScore(ranked, upSet = genes, knownDirection = TRUE,
                     centerScore = TRUE)
  data.frame(
    sample_id = rownames(res),
    score = res$TotalScore,
    stringsAsFactors = FALSE
  )
}
ref_rows <- rbind(
  cbind(set_name = "responder_up",   score_one(up_sig)),
  cbind(set_name = "responder_down", score_one(down_sig)),
  cbind(set_name = "scattered",      score_one(scattered))
)
ref_rows <- ref_rows[order(ref_rows$set_name, ref_rows$sample_id), ]
ref_rows$score <- formatC(ref_rows$score, digits = 17, format = "g")
out_path <- file.path(fixture_dir, "singscore_reference.tsv")
write.table(ref_rows, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "with", nrow(ref_rows), "rows\n")
print(ref_rows)
