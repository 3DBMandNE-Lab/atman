# Reference fit of DEqMS (Bioconductor) on the exact atman CPTAC test
# fixture. Writes deqms_cptac_reference.tsv next to this script so the
# atman integration test can diff against real DEqMS output on identical
# input.
#
# Workflow (mirrors the DEqMS vignette):
#   peptideRaw  -> log2 -> per-sample median-center
#               -> protein-level median summary
#               -> limma (lmFit + eBayes, trend = TRUE)
#               -> DEqMS::spectraCounteBayes(fit, coef_col, count)
#   where `count` is the peptide-count per protein.
#
# Run: Rscript deqms_cptac_reference.R

suppressPackageStartupMessages({
  library(limma)
  library(DEqMS)
  library(matrixStats)
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

# log2 + per-sample median-center. center.median subtracts each sample's
# median log2-intensity (over peptides observed in that sample).
log2_peptides <- log2(intens)
sample_medians <- apply(log2_peptides, 2, median, na.rm = TRUE)
log2_peptides <- sweep(log2_peptides, 2, sample_medians, FUN = "-")
rownames(log2_peptides) <- paste0("P", sprintf("%04d", seq_len(nrow(dat))),
                                  "_", dat$Sequence)

# Protein-level median summary (per sample median across peptides of the
# protein). NA → median() returns NA when all peptides are missing for
# that protein/sample pair.
proteins <- dat$LeadingRazorProtein
uniq <- sort(unique(proteins))
prot_mat <- matrix(
  NA_real_, nrow = length(uniq), ncol = ncol(log2_peptides),
  dimnames = list(uniq, sample_cols)
)
peptide_counts <- integer(length(uniq))
names(peptide_counts) <- uniq
for (p in uniq) {
  rows <- which(proteins == p)
  peptide_counts[p] <- length(rows)
  sub <- log2_peptides[rows, , drop = FALSE]
  if (length(rows) == 1) {
    prot_mat[p, ] <- sub[1, ]
  } else {
    prot_mat[p, ] <- colMedians(sub, na.rm = TRUE)
  }
}

# Design + contrast: ~0 + condition with conditionB - conditionA.
condition <- factor(substr(sample_cols, 1, 1), levels = c("A", "B"))
design <- model.matrix(~0 + condition)
colnames(design) <- c("A", "B")
contrast <- makeContrasts("B-A", levels = design)

fit <- lmFit(prot_mat, design)
fit2 <- contrasts.fit(fit, contrast)
fit3 <- eBayes(fit2, trend = TRUE)

# DEqMS: attach peptide count to the fit and call spectraCounteBayes.
# DEqMS expects `fit$count` on the same order as `fit$coefficients`.
fit3$count <- peptide_counts[rownames(fit3$coefficients)]
fit4 <- spectraCounteBayes(fit3)
tbl <- outputResult(fit4, coef_col = 1)

out <- data.frame(
  assay_id = rownames(tbl),
  logFC_B_minus_A = tbl$logFC,
  se = tbl$logFC / tbl$t,
  t_stat = tbl$t,
  df_total = fit4$df.total[match(rownames(tbl), rownames(fit4$coefficients))],
  p_value = tbl$P.Value,
  adj_p_value = tbl$adj.P.Val,
  sca_t = tbl$sca.t,
  sca_p = tbl$sca.P.Value,
  sca_adj_p = tbl$sca.adj.pval,
  peptide_count = tbl$count,
  stringsAsFactors = FALSE
)
out$class <- ifelse(grepl("_HUMAN_UPS", out$assay_id), "UPS1",
             ifelse(grepl("_YEAST", out$assay_id), "YEAST", "OTHER"))
out <- out[order(out$assay_id), ]

out_path <- file.path(fixture_dir, "deqms_cptac_reference.tsv")
write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "with", nrow(out), "rows\n")

cat("\n=== DEqMS reference summary ===\n")
ups <- out[out$class == "UPS1" & !is.na(out$logFC_B_minus_A), "logFC_B_minus_A"]
yst <- out[out$class == "YEAST" & !is.na(out$logFC_B_minus_A), "logFC_B_minus_A"]
cat(sprintf("UPS1  n=%d  median_logFC_B-A=%.3f  range=[%.3f, %.3f]\n",
            length(ups), median(ups), min(ups), max(ups)))
cat(sprintf("YEAST n=%d  median_logFC_B-A=%.3f  range=[%.3f, %.3f]\n",
            length(yst), median(yst), min(yst), max(yst)))
cat(sprintf("Expected UPS1 logFC(B-A) = %.3f\n", log2(0.74 / 0.25)))
