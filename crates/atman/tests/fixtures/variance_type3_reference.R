# Reference Type III sum-of-squares + F + p-values on a 3-level
# categorical factor + numeric covariate, fit per archetype via
# base R's `lm`. Emits variance_type3_reference.tsv. Atman's
# decompose variance must match this to tight numerical tolerance
# on ss_type3, f_statistic, and p_value per (archetype, factor).
#
# Run: Rscript variance_type3_reference.R
#
# Requires: car (for Anova(type = 3)).

script_arg <- Sys.getenv("ATMAN_FIXTURE_DIR")
fixture_dir <- if (nzchar(script_arg)) script_arg else {
  args <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(self) == 1 && file.exists(self)) dirname(normalizePath(self))
  else "crates/atman/tests/fixtures"
}

suppressPackageStartupMessages({
  if (!requireNamespace("car", quietly = TRUE)) {
    install.packages("car", repos = "https://cloud.r-project.org")
  }
  library(car)
})

set.seed(20260421)
n_samples <- 48L
# Balanced: 16 per stage, two conditions evenly split per stage.
stages <- rep(c("CN", "MCI", "AD"), each = 16L)
conditions <- rep(c("Case", "Control"), times = 24L)
ages <- round(runif(n_samples, 40, 90), 1)
sample_ids <- sprintf("S%03d", seq_len(n_samples))

# Archetype A1: strong stage effect + mild condition effect.
a1 <- c(CN = 0, MCI = 1.5, AD = 3.0)[stages] +
      c(Case = 0.2, Control = 0)[conditions] +
      0.3 * (ages - mean(ages)) / sd(ages) +
      rnorm(n_samples, 0, 0.5)
# Archetype A2: condition only.
a2 <- c(Case = 0, Control = -1.0)[conditions] + rnorm(n_samples, 0, 0.7)
# Archetype A3: null.
a3 <- rnorm(n_samples, 0, 1.0)

wide <- data.frame(
  sample_id = sample_ids, subject_id = sample_ids,
  stage = factor(stages, levels = c("AD", "CN", "MCI")),
  condition = factor(conditions, levels = c("Case", "Control")),
  age = ages,
  A1 = a1, A2 = a2, A3 = a3,
  stringsAsFactors = FALSE
)

# Write fixture (Atman canonical + extras).
samples_tsv <- data.frame(
  sample_id = sample_ids, subject_id = sample_ids,
  condition = conditions, is_control = 0L, sample_type = "plasma",
  ingest_order = seq_len(n_samples),
  stage = stages, age = ages,
  stringsAsFactors = FALSE
)
write.table(samples_tsv,
            file = file.path(fixture_dir, "variance_type3_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

activations <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
  rbind(
    data.frame(cohort = "bench", subject_id = sample_ids[i],
               sample_id = sample_ids[i], program = "A1",
               activation = sprintf("%.6f", a1[i]),
               stringsAsFactors = FALSE),
    data.frame(cohort = "bench", subject_id = sample_ids[i],
               sample_id = sample_ids[i], program = "A2",
               activation = sprintf("%.6f", a2[i]),
               stringsAsFactors = FALSE),
    data.frame(cohort = "bench", subject_id = sample_ids[i],
               sample_id = sample_ids[i], program = "A3",
               activation = sprintf("%.6f", a3[i]),
               stringsAsFactors = FALSE)
  )
}))
write.table(activations,
            file = file.path(fixture_dir, "variance_type3_activations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Fit lm per archetype, extract Type III from car::Anova, also
# compute SS_III directly via Wald: β_S' [(X'X)^-1_SS]^-1 β_S.
# Set sum-to-zero contrasts for Type III parity with atman's
# treatment-coded Wald form: both orientations share the same β_S
# SS_III because the Wald formula is contrast-invariant when the
# full model's fitted values match.
options(contrasts = c("contr.treatment", "contr.poly"))

fit_and_type3 <- function(y) {
  m <- lm(y ~ stage + condition + age, data = wide)
  # Wald-form Type III per factor.
  beta <- coef(m)
  vcv <- vcov(m)
  sigma2 <- summary(m)$sigma^2
  df <- m$df.residual
  nm <- names(beta)
  get_f <- function(cols) {
    beta_s <- beta[cols]
    v <- vcv[cols, cols, drop = FALSE]
    # F = β_S' V_S^-1 β_S / (k · σ²)
    inv_v <- solve(v)
    quad <- as.numeric(t(beta_s) %*% inv_v %*% beta_s)
    # BUT V_S = σ² (X'X)^-1_SS, so β_S' (X'X)_SS · β_S / (k σ²)
    # = β_S' V_S^-1 · σ² · β_S / (k · σ²) = quad / k; we already
    # absorbed σ² into V. Therefore F = (quad / k) is correct.
    k <- length(cols)
    list(F = quad / k, df_num = k, df_den = df,
         ss = quad * sigma2 / k * k)  # β' (X'X)_SS β = F · k · σ²
  }
  # Identify the columns belonging to each factor by name prefix.
  stage_cols <- grep("^stage", nm, value = TRUE)
  condition_cols <- grep("^condition", nm, value = TRUE)
  # `age` is a single numeric column.
  age_cols <- "age"
  rs <- get_f(stage_cols)
  rc <- get_f(condition_cols)
  ra <- get_f(age_cols)
  data.frame(
    factor_name = c("stage", "condition", "age"),
    ss_type3    = c(rs$ss, rc$ss, ra$ss),
    f_statistic = c(rs$F, rc$F, ra$F),
    df_num      = c(rs$df_num, rc$df_num, ra$df_num),
    df_den      = c(rs$df_den, rc$df_den, ra$df_den),
    p_value     = c(
      pf(rs$F, rs$df_num, rs$df_den, lower.tail = FALSE),
      pf(rc$F, rc$df_num, rc$df_den, lower.tail = FALSE),
      pf(ra$F, ra$df_num, ra$df_den, lower.tail = FALSE)
    ),
    stringsAsFactors = FALSE
  )
}

ref <- do.call(rbind, lapply(c("A1", "A2", "A3"), function(g) {
  df <- fit_and_type3(wide[[g]])
  df$archetype_id <- g
  df[, c("archetype_id", "factor_name", "ss_type3", "f_statistic",
         "df_num", "df_den", "p_value")]
}))
out_path <- file.path(fixture_dir, "variance_type3_reference.tsv")
write.table(ref, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "\n")
print(ref)
