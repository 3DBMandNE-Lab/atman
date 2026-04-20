# Reference Sidak-adjusted pairwise contrast p-values on a 3-level
# categorical factor fit via OLS. Emits posthoc_sidak_reference.tsv
# with one row per protein × contrast. Atman's de --post-hoc sidak
# must match this TSV's posthoc_p and posthoc_adj_p within tight
# numerical tolerance.
#
# Factor: stage ∈ {CN, MCI, AD}  (alphabetical reference: AD).
# Contrasts: MCI-CN, AD-CN, AD-MCI.
# Design: y ~ stage + age
# Sidak adjustment: p_adj = 1 - (1 - p)^m, m = 3.
#
# Run: Rscript posthoc_sidak_reference.R

script_arg <- Sys.getenv("ATMAN_FIXTURE_DIR")
fixture_dir <- if (nzchar(script_arg)) script_arg else {
  args <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(self) == 1 && file.exists(self)) dirname(normalizePath(self))
  else "crates/atman/tests/fixtures"
}

set.seed(20260420)
n_samples <- 36L
stages <- rep(c("CN", "MCI", "AD"), each = 12L)
ages <- round(runif(n_samples, 50, 80), 1)
sample_ids <- sprintf("S%03d", seq_len(n_samples))

# Two proteins: responder has a clear stage effect; null has none.
# Fixed deterministic noise via R's Gaussian RNG (seed above).
responder_effect <- c(CN = 0.0, MCI = 1.2, AD = 2.1)
responder <- responder_effect[stages] + 0.4 * (ages - mean(ages)) / sd(ages) +
             rnorm(n_samples, 0, 0.5)
null_protein <- rnorm(n_samples, 0, 1.0)

wide <- data.frame(
  sample_id = sample_ids,
  stage = factor(stages, levels = c("AD", "CN", "MCI")),
  age = ages,
  RESPONDER = responder,
  NULLP = null_protein,
  stringsAsFactors = FALSE
)
# Emit per-sample metadata for atman (canonical sample sheet layout +
# stage column appended).
samples_tsv <- data.frame(
  sample_id = sample_ids,
  subject_id = sample_ids,
  condition = "Case",
  is_control = 0L,
  sample_type = "plasma",
  ingest_order = seq_len(n_samples),
  stage = stages,
  age = ages,
  stringsAsFactors = FALSE
)
write.table(samples_tsv, file = file.path(fixture_dir, "posthoc_sidak_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Proteins catalog.
prot <- data.frame(
  platform = "olink_explore_ngs",
  assay_id = c("R001", "N001"),
  uniprot = c("Q00001", "Q00002"),
  gene_symbol = c("RESPONDER", "NULLP"),
  panel = "P1",
  panel_lot = "",
  stringsAsFactors = FALSE
)
write.table(prot, file = file.path(fixture_dir, "posthoc_sidak_proteins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# qc_measurements.tsv long format.
qc <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
  sid <- sample_ids[i]
  rbind(
    data.frame(platform = "olink_explore_ngs", sample_id = sid, assay_id = "R001",
               gene_symbol = "RESPONDER", panel = "P1",
               npx_source_str = sprintf("%.6f", responder[i]),
               abundance = sprintf("%.6f", responder[i]),
               abundance_raw = sprintf("%.6f", responder[i]),
               abundance_unit = "log2_npx",
               qc_sample = "PASS", qc_assay = "PASS",
               detection_limit = "", below_lod = 0L, dropped_by_qc = 0L,
               plate_id = "", panel_lot = "",
               ingest_order = 2 * (i - 1) + 1,
               stringsAsFactors = FALSE),
    data.frame(platform = "olink_explore_ngs", sample_id = sid, assay_id = "N001",
               gene_symbol = "NULLP", panel = "P1",
               npx_source_str = sprintf("%.6f", null_protein[i]),
               abundance = sprintf("%.6f", null_protein[i]),
               abundance_raw = sprintf("%.6f", null_protein[i]),
               abundance_unit = "log2_npx",
               qc_sample = "PASS", qc_assay = "PASS",
               detection_limit = "", below_lod = 0L, dropped_by_qc = 0L,
               plate_id = "", panel_lot = "",
               ingest_order = 2 * i,
               stringsAsFactors = FALSE)
  )
}))
write.table(qc, file = file.path(fixture_dir, "posthoc_sidak_qc_measurements.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, file = file.path(fixture_dir, "posthoc_sidak_measurements.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# R reference fits: lm per protein, pairwise contrasts via
# coefficient combinations. Sidak adjust across 3 contrasts.
fit_and_contrast <- function(y) {
  m <- lm(y ~ stage + age, data = wide)
  co <- coef(m)
  vcv <- vcov(m)
  df <- m$df.residual
  # Under R's default "AD" reference (alphabetically first level
  # given levels=c("AD","CN","MCI")):
  #   β_CN  = coef "stageCN"   -> mean(CN) - mean(AD)
  #   β_MCI = coef "stageMCI"  -> mean(MCI) - mean(AD)
  # Atman's contrast orientation is `level_a - level_b` in the
  # label "MCI-CN" → estimate = mean(MCI) - mean(CN) = β_MCI - β_CN.
  eval_contrast <- function(c) {
    est <- sum(c * co)
    var_est <- as.numeric(t(c) %*% vcv %*% c)
    se <- sqrt(var_est)
    t <- est / se
    p <- 2 * pt(-abs(t), df)
    list(est = est, se = se, t = t, p = p, df = df)
  }
  c_names <- names(co)
  zero <- setNames(numeric(length(c_names)), c_names)
  build <- function(a, b) {
    c <- zero
    # Atman convention: estimate = a - b, so weight +1 at level_a
    # column and -1 at level_b column (with reference level
    # contributing 0).
    if (a == "AD" && b != "AD") c[[paste0("stage", b)]] <- -1
    else if (b == "AD" && a != "AD") c[[paste0("stage", a)]] <-  1
    else if (a != "AD" && b != "AD") {
      c[[paste0("stage", a)]] <-  1
      c[[paste0("stage", b)]] <- -1
    }
    c
  }
  list(
    `MCI-CN`  = eval_contrast(build("MCI", "CN")),
    `AD-CN`   = eval_contrast(build("AD", "CN")),
    `AD-MCI`  = eval_contrast(build("AD", "MCI"))
  )
}

m_contrasts <- 3L
sidak_adj <- function(p) 1 - (1 - pmin(pmax(p, 0), 1))^m_contrasts

rows <- list()
for (gene in c("RESPONDER", "NULLP")) {
  y <- wide[[gene]]
  cr <- fit_and_contrast(y)
  for (lbl in names(cr)) {
    r <- cr[[lbl]]
    rows[[length(rows) + 1L]] <- data.frame(
      gene_symbol = gene,
      comparison = lbl,
      estimate = r$est,
      se = r$se,
      t = r$t,
      df = r$df,
      posthoc_p = r$p,
      posthoc_adj_p = sidak_adj(r$p),
      stringsAsFactors = FALSE
    )
  }
}
ref <- do.call(rbind, rows)
out_path <- file.path(fixture_dir, "posthoc_sidak_reference.tsv")
write.table(ref, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "with", nrow(ref), "rows\n")
print(ref)
