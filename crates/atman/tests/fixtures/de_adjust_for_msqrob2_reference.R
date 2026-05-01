#!/usr/bin/env Rscript
# crates/atman/tests/fixtures/de_adjust_for_msqrob2_reference.R
#
# Builds a peptide-level QFeatures object from the de_adjust_for_input
# fixture (protein-level data replicated into 3 synthetic peptides per
# protein with fixed offsets -0.1 / 0.0 / +0.1), fits msqrob2 with
# formula ~ condition + nmf_program_01 + nmf_program_02, and emits a
# reference TSV for the K-phase parity test.
#
# Peptide offsets are deterministic: the same offsets are hard-coded in the
# atman peptide fixture files so both sides see identical input data.
#
# Pinned: msqrob2 1.16.0, R 4.5.0
# Coefficient tested: "conditionb" (condition b vs a; sign = mean_b - mean_a)
#
# Usage:
#   Rscript de_adjust_for_msqrob2_reference.R \
#     <input_dir> <covariates_tsv> <output_tsv>

suppressPackageStartupMessages({
  library(msqrob2)
  library(QFeatures)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript de_adjust_for_msqrob2_reference.R <input_dir> <covariates_tsv> <output_tsv>")
}
input_dir      <- args[1]
covariates_tsv <- args[2]
output_tsv     <- args[3]

# ── helpers ───────────────────────────────────────────────────────────────────
read_tsv_skip_comment <- function(path) {
  raw  <- readLines(path)
  body <- raw[!startsWith(raw, "#")]
  read.table(text = paste(body, collapse = "\n"),
             sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "")
}

# ── read canonical inputs ─────────────────────────────────────────────────────
samples      <- read_tsv_skip_comment(file.path(input_dir, "samples.tsv"))
measurements <- read_tsv_skip_comment(file.path(input_dir, "measurements.tsv"))
covariates   <- read_tsv_skip_comment(covariates_tsv)

# Non-control samples only.
samples <- samples[samples$is_control == 0 & !is.na(samples$condition), ]
sample_ids <- samples$sample_id  # 12 samples in the fixture

# Protein-level wide matrix (genes × samples).
genes   <- sort(unique(measurements$gene_symbol))
n_genes <- length(genes)
n_samp  <- length(sample_ids)

prot_mat <- matrix(NA_real_, nrow = n_genes, ncol = n_samp,
                   dimnames = list(genes, sample_ids))
for (i in seq_len(nrow(measurements))) {
  g   <- measurements$gene_symbol[i]
  sid <- measurements$sample_id[i]
  if (sid %in% sample_ids && g %in% genes) {
    prot_mat[g, sid] <- as.numeric(measurements$abundance[i])
  }
}

# ── synthetic peptide-level matrix ────────────────────────────────────────────
# 3 peptides per protein: offsets -0.1, 0.0, +0.1 (deterministic; matches
# atman fixture files in de_adjust_for_input/).
pep_offsets <- c(-0.1, 0.0, 0.1)
n_pep       <- length(pep_offsets)
pep_names   <- paste0(rep(genes, each = n_pep),
                      "_pep", rep(seq_len(n_pep), times = n_genes))
pep_protein <- rep(genes, each = n_pep)

pep_mat <- matrix(NA_real_, nrow = n_genes * n_pep, ncol = n_samp,
                  dimnames = list(pep_names, sample_ids))
for (g in genes) {
  for (k in seq_along(pep_offsets)) {
    row_name <- paste0(g, "_pep", k)
    pep_mat[row_name, ] <- prot_mat[g, ] + pep_offsets[k]
  }
}

# ── build QFeatures ───────────────────────────────────────────────────────────
rownames(covariates) <- covariates$sample_id
cov_ordered <- covariates[sample_ids, , drop = FALSE]

condition <- factor(samples$condition[match(sample_ids, samples$sample_id)],
                    levels = c("a", "b"))

cd <- DataFrame(
  sample_id      = sample_ids,
  condition      = condition,
  nmf_program_01 = cov_ordered$nmf_program_01,
  nmf_program_02 = cov_ordered$nmf_program_02,
  row.names      = sample_ids
)

rd_pep <- DataFrame(
  peptide_id = pep_names,
  Proteins   = pep_protein,
  row.names  = pep_names
)

se_pep <- SummarizedExperiment(
  assays  = list(peptideIntensity = pep_mat),
  rowData = rd_pep,
  colData = cd
)

pe <- QFeatures(list(peptideIntensity = se_pep), colData = cd)

# ── aggregate peptides → protein level ───────────────────────────────────────
# Use robust summarisation (default in msqrob2 vignette). With 3 identical-
# offset peptides the robust summary coincides with the mean, which equals
# the original protein abundance.
pe <- aggregateFeatures(
  pe, i = "peptideIntensity", fcol = "Proteins", name = "protein",
  fun = MsCoreUtils::robustSummary, na.rm = TRUE
)

# ── msqrob fit with covariates ────────────────────────────────────────────────
pe <- msqrob(
  object  = pe,
  i       = "protein",
  formula = ~ condition + nmf_program_01 + nmf_program_02
)

# ── hypothesis test: conditionb ───────────────────────────────────────────────
L  <- makeContrast("conditionb = 0", parameterNames = "conditionb")
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

# ── extract results ───────────────────────────────────────────────────────────
rd_out <- rowData(pe[["protein"]])
tbl    <- rd_out[["conditionb"]]
stopifnot(!is.null(tbl))

# Gene symbols from row names of the protein assay.
protein_ids <- rownames(pe[["protein"]])
# protein_ids are gene symbols (they were used as row keys in prot_mat).
gene_symbols <- protein_ids

out <- data.frame(
  gene_symbol = gene_symbols,
  log_fc      = tbl$logFC,
  t_stat      = tbl$t,
  p_value     = tbl$pval,
  q_value     = tbl$adjPval,
  stringsAsFactors = FALSE
)

# Sort deterministically by gene_symbol to match atman's BTreeMap ordering.
out <- out[order(out$gene_symbol), ]

# ── emit reference TSV ────────────────────────────────────────────────────────
msqrob2_ver <- as.character(packageVersion("msqrob2"))
r_ver       <- paste(R.version$major, R.version$minor, sep = ".")
header_comment <- paste0(
  "# generated by de_adjust_for_msqrob2_reference.R",
  " msqrob2=", msqrob2_ver,
  " R=", r_ver,
  " coef=conditionb sign=mean_b_minus_mean_a"
)

con <- file(output_tsv, open = "w")
writeLines(header_comment, con)
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

cat(sprintf("wrote %s (%d genes)\n", output_tsv, nrow(out)))
cat(sprintf("msqrob2 %s, R %s\n", msqrob2_ver, r_ver))
cat(sprintf("formula: ~ condition + nmf_program_01 + nmf_program_02\n"))
cat(sprintf("coef tested: conditionb (mean_b - mean_a)\n"))
