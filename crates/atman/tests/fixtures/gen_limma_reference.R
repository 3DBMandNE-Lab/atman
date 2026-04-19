#!/usr/bin/env Rscript
# Generate limma reference TSVs for atman's fixture-based tests.
#
# Run manually when limma changes (not executed by CI):
#   cd crates/atman/tests/fixtures
#   Rscript gen_limma_reference.R
#
# Inputs:  limma_fixture_samples.tsv, limma_fixture_proteins.tsv,
#          limma_fixture_measurements.tsv
# Outputs: limma_reference_notrend.tsv, limma_reference_trend.tsv
#
# ---------------------------------------------------------------------------
# SIGN CONVENTION NOTE (load-bearing for Task 14)
# ---------------------------------------------------------------------------
# atman's `A-B` comparison semantics report `mean_a - mean_b`. This R
# reference uses limma's `~groups` design with `A` as the reference level,
# so the `groupB` coefficient (`fit$coef[, 2]`) equals `mean_b - mean_a` --
# the opposite sign from atman's convention.
#
# The reference TSVs below are written VERBATIM from limma (positive
# `mean_b - mean_a` for `t` and the corresponding `P_Value`). Task 12's
# Rust code applies contrast `[0, -1]` so `atman de --test limma
# --groups A-B` emits `mean_a - mean_b` (sign-flipped relative to these
# fixtures).
#
# Task 14's test code is responsible for reconciling the sign: either
# compare by absolute value, or negate R's `t` at read time. Do NOT flip
# at fixture-generation time -- the fixture files should be verbatim R
# limma output for auditability.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(limma))

set.seed(20260418)

# --- Synthetic dataset generator -------------------------------------------
# 100 features x 20 samples, two-group design (A vs B, 10 each). Effects are
# planted on rows 1..50, null on 51..100.
n_features <- 100
n_samples  <- 20
group_sizes <- c(10, 10)
groups <- factor(c(rep("A", group_sizes[1]), rep("B", group_sizes[2])),
                 levels = c("A", "B"))

true_effect <- c(rep(1.0, 50), rep(0.0, 50))
# Per-feature residual variance drawn from a scaled chi-square / df_prior = 6
df_prior_true <- 6
s2_prior_true <- 0.25
feat_sigma2 <- s2_prior_true * df_prior_true / rchisq(n_features, df_prior_true)
feat_sigma  <- sqrt(feat_sigma2)

y <- matrix(0, nrow = n_features, ncol = n_samples)
rownames(y) <- sprintf("F%03d", seq_len(n_features))
colnames(y) <- sprintf("S%02d", seq_len(n_samples))
# Shift the per-feature baseline to a positive range so log(Amean) is
# well-defined for the parametric quadratic trend fit. 5.0 matches the
# typical log2-NPX abundance scale and keeps Amean strictly positive.
baseline <- 5.0
for (i in seq_len(n_features)) {
    shift_b <- true_effect[i]
    x_a <- rnorm(group_sizes[1], mean = baseline,           sd = feat_sigma[i])
    x_b <- rnorm(group_sizes[2], mean = baseline + shift_b, sd = feat_sigma[i])
    y[i, ] <- c(x_a, x_b)
}

# --- Write atman canonical TSVs for the atman binary to consume ------------
fixture_dir <- "."

samples <- data.frame(
    sample_id    = colnames(y),
    subject_id   = colnames(y),
    condition    = as.character(groups),
    is_control   = 0,
    sample_type  = "plasma",
    ingest_order = seq_len(n_samples),
    stringsAsFactors = FALSE
)
write.table(samples,
            file.path(fixture_dir, "limma_fixture_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

proteins <- data.frame(
    platform    = "olink_explore_ngs",
    assay_id    = rownames(y),
    uniprot     = sprintf("Q%05d", seq_len(n_features)),
    gene_symbol = sprintf("GENE%03d", seq_len(n_features)),
    panel       = "P1",
    panel_lot   = "",
    stringsAsFactors = FALSE
)
write.table(proteins,
            file.path(fixture_dir, "limma_fixture_proteins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

meas_rows <- list()
ingest_order <- 0
for (i in seq_len(n_features)) {
    for (j in seq_len(n_samples)) {
        ingest_order <- ingest_order + 1
        v <- y[i, j]
        meas_rows[[ingest_order]] <- sprintf(
            paste(c(
                "olink_explore_ngs", colnames(y)[j], rownames(y)[i],
                sprintf("GENE%03d", i), "P1", sprintf("%.6f", v),
                sprintf("%.6f", v), sprintf("%.6f", v), "log2_npx",
                "PASS", "PASS", "", "0", "0", "", "", as.character(ingest_order)
            ), collapse = "\t")
        )
    }
}
meas_text <- paste0(
    paste(c("platform","sample_id","assay_id","gene_symbol","panel",
            "npx_source_str","abundance","abundance_raw","abundance_unit",
            "qc_sample","qc_assay","detection_limit","below_lod",
            "dropped_by_qc","plate_id","panel_lot","ingest_order"), collapse = "\t"),
    "\n",
    paste(unlist(meas_rows), collapse = "\n"), "\n")
writeLines(meas_text, file.path(fixture_dir, "limma_fixture_measurements.tsv"))

# --- Fit limma --------------------------------------------------------------
design <- model.matrix(~ groups)
fit <- lmFit(y, design)

# Reference 1: no trend.
fit_notrend <- eBayes(fit, trend = FALSE, robust = TRUE)
ref_notrend <- data.frame(
    feature      = rownames(y),
    s2_trend     = NA_real_,
    s2_prior     = fit_notrend$s2.prior,
    s2_posterior = fit_notrend$s2.post,
    df_prior     = fit_notrend$df.prior,
    df_total     = fit_notrend$df.total,
    t            = fit_notrend$t[, 2],
    P_Value      = fit_notrend$p.value[, 2],
    F_statistic  = fit_notrend$F,
    F_p_value    = fit_notrend$F.p.value,
    stringsAsFactors = FALSE
)
write.table(ref_notrend,
            file.path(fixture_dir, "limma_reference_notrend.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Reference 2: trend, passing a precomputed parametric quadratic so the Rust
# port's trend matches R's byte-for-byte.
Amean <- rowMeans(y)
lm_trend <- lm(log(fit$sigma^2) ~ log(Amean) + I(log(Amean)^2))
s2_trend <- exp(predict(lm_trend))
fit_trend <- eBayes(fit, trend = s2_trend, robust = TRUE)
ref_trend <- data.frame(
    feature      = rownames(y),
    s2_trend     = s2_trend,
    s2_prior     = fit_trend$s2.prior,
    s2_posterior = fit_trend$s2.post,
    df_prior     = fit_trend$df.prior,
    df_total     = fit_trend$df.total,
    t            = fit_trend$t[, 2],
    P_Value      = fit_trend$p.value[, 2],
    F_statistic  = fit_trend$F,
    F_p_value    = fit_trend$F.p.value,
    stringsAsFactors = FALSE
)
write.table(ref_trend,
            file.path(fixture_dir, "limma_reference_trend.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote limma_reference_{notrend,trend}.tsv and canonical fixture TSVs.\n")
cat("sessionInfo():\n")
print(sessionInfo())
