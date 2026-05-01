#!/usr/bin/env Rscript
# crates/atman/tests/fixtures/de_adjust_for_msqrob2_reference.R
#
# Runs msqrob2 on peptide-level data matching the atman fixture
# (de_adjust_for_input: 12 samples × 30 proteins × 3 peptides per protein,
# with fixed offsets -0.1/0.0/+0.1 from the protein abundance).
#
# Uses msqrobAggregate with:
#   formula  = ~ condition + nmf_program_01 + nmf_program_02 + (1|peptide_id)
#   fcol     = "Proteins"   (grouping variable → protein-level models)
#
# This dispatches msqrobLmer per protein on the stacked peptide-level
# observations, matching atman's LMM with random peptide intercept.
# The reported coefficient is conditionb = mean_b − mean_a.
#
# Pinned: msqrob2 1.16.0, R 4.5.0
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
samples       <- read_tsv_skip_comment(file.path(input_dir, "samples.tsv"))
measurements  <- read_tsv_skip_comment(file.path(input_dir, "measurements.tsv"))
covariates    <- read_tsv_skip_comment(covariates_tsv)
peptides_meta <- read.table(file.path(input_dir, "peptides.tsv"),
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pep_meas      <- read.table(file.path(input_dir, "peptide_measurements.tsv"),
                            sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Non-control samples only.
samples    <- samples[samples$is_control == 0 & !is.na(samples$condition), ]
sample_ids <- samples$sample_id  # 12 samples

# ── peptide-level intensity matrix ────────────────────────────────────────────
pep_ids <- sort(unique(pep_meas$peptide_id))
n_pep   <- length(pep_ids)
n_samp  <- length(sample_ids)

pep_mat <- matrix(NA_real_, nrow = n_pep, ncol = n_samp,
                  dimnames = list(pep_ids, sample_ids))
for (i in seq_len(nrow(pep_meas))) {
  pid <- pep_meas$peptide_id[i]
  sid <- pep_meas$sample_id[i]
  if (pid %in% pep_ids && sid %in% sample_ids) {
    pep_mat[pid, sid] <- as.numeric(pep_meas$abundance[i])
  }
}

# ── colData: sample metadata + covariates ─────────────────────────────────────
rownames(covariates) <- covariates$sample_id
cov_ordered <- covariates[sample_ids, , drop = FALSE]
condition   <- factor(samples$condition[match(sample_ids, samples$sample_id)],
                      levels = c("a", "b"))
cd <- DataFrame(
  sample_id      = sample_ids,
  condition      = condition,
  nmf_program_01 = cov_ordered$nmf_program_01,
  nmf_program_02 = cov_ordered$nmf_program_02,
  row.names      = sample_ids
)

# ── rowData: peptide_id (random effect) and Proteins (grouping) ───────────────
# Map: peptide_id → assay_id → gene_symbol.
pid_to_aid  <- setNames(peptides_meta$assay_id, peptides_meta$peptide_id)
prot_meta   <- unique(measurements[, c("assay_id", "gene_symbol")])
aid_to_gene <- setNames(prot_meta$gene_symbol, prot_meta$assay_id)
pep_gene    <- aid_to_gene[pid_to_aid[pep_ids]]

rd_pep <- DataFrame(
  peptide_id = pep_ids,    # random-effect grouping variable (in formula)
  Proteins   = pep_gene,   # protein-level grouping variable (fcol)
  row.names  = pep_ids
)

# ── QFeatures object ──────────────────────────────────────────────────────────
se_pep <- SummarizedExperiment(
  assays  = list(peptideIntensity = pep_mat),
  rowData = rd_pep,
  colData = cd
)
pe <- QFeatures(list(peptideIntensity = se_pep), colData = cd)

# ── msqrobAggregate: peptide-level LMM, one model per protein ────────────────
# Formula: fixed effects (condition, covariates) + random peptide intercept.
# msqrobAggregate dispatches msqrobLmer per protein (featureGroups = Proteins)
# and stores models in the aggregated assay's rowData.
# robust = FALSE matches atman's --robust false flag in the test.
pe <- msqrobAggregate(
  object    = pe,
  i         = "peptideIntensity",
  formula   = ~ condition + nmf_program_01 + nmf_program_02 + (1|peptide_id),
  fcol      = "Proteins",
  name      = "protein",
  robust    = FALSE
)

# ── hypothesis test: conditionb ───────────────────────────────────────────────
L  <- makeContrast("conditionb = 0", parameterNames = "conditionb")
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)

# ── extract results ───────────────────────────────────────────────────────────
rd_out <- rowData(pe[["protein"]])
tbl    <- rd_out[["conditionb"]]
stopifnot(!is.null(tbl))

protein_ids  <- rownames(pe[["protein"]])  # gene symbols (used as row keys)
out <- data.frame(
  gene_symbol = protein_ids,
  log_fc      = tbl$logFC,
  t_stat      = tbl$t,
  p_value     = tbl$pval,
  q_value     = tbl$adjPval,
  stringsAsFactors = FALSE
)
out <- out[order(out$gene_symbol), ]

# ── emit reference TSV ────────────────────────────────────────────────────────
msqrob2_ver <- as.character(packageVersion("msqrob2"))
r_ver       <- paste(R.version$major, R.version$minor, sep = ".")
header_comment <- paste0(
  "# generated by de_adjust_for_msqrob2_reference.R",
  " msqrob2=", msqrob2_ver,
  " R=", r_ver,
  " coef=conditionb sign=mean_b_minus_mean_a",
  " model=peptide_level_lmm_msqrobAggregate"
)

con <- file(output_tsv, open = "w")
writeLines(header_comment, con)
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

cat(sprintf("wrote %s (%d genes)\n", output_tsv, nrow(out)))
cat(sprintf("msqrob2 %s, R %s\n", msqrob2_ver, r_ver))
cat(sprintf("formula: ~ condition + nmf_program_01 + nmf_program_02 + (1|peptide_id)\n"))
cat(sprintf("coef tested: conditionb (mean_b - mean_a)\n"))
cat(sprintf("logFC range: [%.4f, %.4f]\n",
            min(out$log_fc, na.rm = TRUE), max(out$log_fc, na.rm = TRUE)))
