# Reference fit of msqrob2 on the exact atman CPTAC test fixture.
# Mirrors the canonical msqrob2 vignette workflow (log2 → normalize
# median → robust peptide-to-protein summary → msqrob fit with
# `~condition`). Writes msqrob2_cptac_reference.tsv next to this
# script so the atman integration test can diff against real msqrob2
# output on identical input.
#
# Run: Rscript msqrob2_cptac_reference.R
#
# Requires: msqrob2 (Bioconductor, >= 1.16), QFeatures.

suppressPackageStartupMessages({
  library(msqrob2)
  library(QFeatures)
  library(SummarizedExperiment)
})

script_arg <- Sys.getenv("ATMAN_FIXTURE_DIR")
fixture_dir <- if (nzchar(script_arg)) script_arg else {
  args <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(self) == 1 && file.exists(self)) dirname(normalizePath(self))
  else "crates/atman/tests/fixtures"
}
fixture <- file.path(fixture_dir, "cptac_peptides30_6ab.tsv")
stopifnot(file.exists(fixture))

dat <- read.delim(fixture, stringsAsFactors = FALSE, check.names = FALSE)
sample_cols <- c(paste0("A", 1:9), paste0("B", 1:9))
intens <- as.matrix(dat[, sample_cols])
mode(intens) <- "numeric"
intens[intens == 0] <- NA
rownames(intens) <- paste0("P", sprintf("%04d", seq_len(nrow(dat))), "_", dat$Sequence)

rd <- DataFrame(
  peptide_id = rownames(intens),
  Sequence = dat$Sequence,
  Proteins = dat$LeadingRazorProtein
)
cd <- DataFrame(
  sample_id = sample_cols,
  condition = factor(substr(sample_cols, 1, 1), levels = c("A", "B"))
)
rownames(cd) <- sample_cols

se <- SummarizedExperiment(
  assays = list(peptideRaw = intens),
  rowData = rd,
  colData = cd
)
pe <- QFeatures(list(peptideRaw = se), colData = cd)

# Canonical vignette steps.
pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
pe <- normalize(pe, i = "peptideLog", name = "peptideNorm",
                method = "center.median")

# Require each peptide to be observed in ≥ 1 sample of each condition.
pe_norm <- pe[["peptideNorm"]]
conds <- colData(pe_norm)$condition
a_idx <- which(conds == "A")
b_idx <- which(conds == "B")
a_counts <- rowSums(!is.na(assay(pe_norm)[, a_idx, drop = FALSE]))
b_counts <- rowSums(!is.na(assay(pe_norm)[, b_idx, drop = FALSE]))
keep_pep <- a_counts >= 1 & b_counts >= 1

# Require each protein to retain ≥ 2 such peptides (median-polish
# needs multiple rows).
kept <- pe_norm[keep_pep, ]
prot_counts <- table(rowData(kept)$Proteins)
keep_prot <- names(prot_counts[prot_counts >= 2])
kept <- kept[rowData(kept)$Proteins %in% keep_prot, ]
pe <- addAssay(pe, kept, name = "peptideFiltered")

# Median summary (simpler than robustSummary, avoids the internal
# model.matrix call that fails when a protein has only one condition
# represented even after filtering).
pe <- aggregateFeatures(
  pe, i = "peptideFiltered", fcol = "Proteins", name = "protein",
  fun = matrixStats::colMedians, na.rm = TRUE
)

pe <- msqrob(object = pe, i = "protein", formula = ~condition)
L <- makeContrast("conditionB = 0", parameterNames = "conditionB")
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

rd_out <- rowData(pe[["protein"]])
tbl <- rd_out[["conditionB"]]
stopifnot(!is.null(tbl))

out <- data.frame(
  assay_id = rownames(pe[["protein"]]),
  logFC_B_minus_A = tbl$logFC,
  se = tbl$se,
  t_stat = tbl$t,
  df = tbl$df,
  p_value = tbl$pval,
  adj_p_value = tbl$adjPval,
  stringsAsFactors = FALSE
)
out$class <- ifelse(grepl("_HUMAN_UPS", out$assay_id), "UPS1",
             ifelse(grepl("_YEAST", out$assay_id), "YEAST", "OTHER"))
out <- out[order(out$assay_id), ]

out_path <- file.path(fixture_dir, "msqrob2_cptac_reference.tsv")
write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "with", nrow(out), "rows\n")

cat("\n=== msqrob2 (protein-summary workflow) ===\n")
ups <- out[out$class == "UPS1" & !is.na(out$logFC_B_minus_A), "logFC_B_minus_A"]
yst <- out[out$class == "YEAST" & !is.na(out$logFC_B_minus_A), "logFC_B_minus_A"]
cat(sprintf("UPS1  n=%d  median=%.3f  range=[%.3f, %.3f]\n",
            length(ups), median(ups), min(ups), max(ups)))
cat(sprintf("YEAST n=%d  median=%.3f  range=[%.3f, %.3f]\n",
            length(yst), median(yst), min(yst), max(yst)))
cat(sprintf("Expected UPS1 logFC(B-A) = %.3f\n", log2(0.74 / 0.25)))
