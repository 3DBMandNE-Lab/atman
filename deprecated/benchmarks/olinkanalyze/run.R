#!/usr/bin/env Rscript
# OlinkAnalyze comparator on the Dube dataset.
#
# Runs olink_ttest paired at subject level across four contrasts, writes
# results to benchmarks/out/olinkanalyze_de.tsv in the shared schema.
# Also runs olink_lmer as a secondary repeated-measures check.

suppressPackageStartupMessages({
    library(OlinkAnalyze)
    library(dplyr)
    library(readr)
    library(tidyr)
})

# Run from repo root (set via bash wrapper). Guard in case of direct invocation.
if (!file.exists("Cargo.toml")) {
    stop("run from repo root; `Cargo.toml` not found in getwd() = ", getwd())
}

out_dir <- "benchmarks/out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load Dube NPX via OlinkAnalyze's own parser --------------------------
npx_files <- c(
    "example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv",
    "example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv"
)

load_npx <- function(path) {
    raw <- read_NPX(path)
    raw
}
npx_raw <- bind_rows(lapply(npx_files, load_npx))

# Parse Dube SampleID -> participant + condition.
parse_dube <- function(sid) {
    m <- regmatches(sid, regexec("^SSNA-([0-9]+)-(PR1|PT1|PR2|PT2)$", sid))
    out <- do.call(rbind, lapply(m, function(x) if (length(x) == 3) x[2:3] else c(NA, NA)))
    colnames(out) <- c("participant", "condition")
    as.data.frame(out, stringsAsFactors = FALSE)
}
ids <- parse_dube(npx_raw$SampleID)
npx <- npx_raw %>%
    mutate(participant = ids$participant, condition = ids$condition) %>%
    filter(!is.na(condition))           # drop controls
# Drop assay warnings per Dube QC rule.
npx <- npx %>% filter(Assay_Warning == "PASS", QC_Warning == "PASS")

# ---- Paired t-tests across the four contrasts -----------------------------
contrasts <- list(
    "PT1-PR1" = c("PT1", "PR1"),
    "PT2-PR2" = c("PT2", "PR2"),
    "PT2-PT1" = c("PT2", "PT1"),
    "PR2-PR1" = c("PR2", "PR1")
)

t_start <- Sys.time()
all_de <- list()
for (cname in names(contrasts)) {
    pair <- contrasts[[cname]]
    sub <- npx %>% filter(condition %in% pair) %>%
        mutate(condition = factor(condition, levels = pair))
    # olink_ttest is paired when the variable is a factor with 2 levels and
    # the data carries repeated SubjectID — rename participant -> SubjectID.
    sub <- sub %>% rename(SubjectID = participant)
    res <- olink_ttest(df = sub, variable = "condition",
                       pair_id = "SubjectID", alternative = "two.sided")
    # Satterthwaite by default; BH within-contrast follows.
    res <- res %>% mutate(contrast = cname)
    all_de[[cname]] <- res
}
t_end <- Sys.time()
runtime_s <- as.numeric(difftime(t_end, t_start, units = "secs"))

de <- bind_rows(all_de) %>%
    transmute(
        contrast  = contrast,
        gene      = Assay,
        mean_diff = estimate,
        pvalue    = p.value,
        qvalue    = Adjusted_pval,
        tool      = "OlinkAnalyze_ttest",
        runtime_s = runtime_s,
        peak_rss_mb = NA_real_
    )
write_tsv(de, file.path(out_dir, "olinkanalyze_de.tsv"))
cat(sprintf("olinkanalyze_ttest DE rows: %d  runtime=%.2fs\n", nrow(de), runtime_s))

# ---- Secondary: olink_lmer (repeated-measures LMM) ------------------------
# Guard — olink_lmer requires >=2 levels of the variable; we fit across all
# four conditions and extract the acute vs rest contrast at post-hoc stage.
t_start <- Sys.time()
npx_lmer <- npx %>% rename(SubjectID = participant) %>%
    mutate(condition = factor(condition, levels = c("PR1", "PT1", "PR2", "PT2")))
lmer_res <- tryCatch({
    olink_lmer(df = npx_lmer, variable = "condition", random = "SubjectID")
}, error = function(e) { message("olink_lmer failed: ", e$message); NULL })
t_end <- Sys.time()
lmer_runtime <- as.numeric(difftime(t_end, t_start, units = "secs"))
if (!is.null(lmer_res)) {
    write_tsv(lmer_res, file.path(out_dir, "olinkanalyze_lmer.tsv"))
    cat(sprintf("olinkanalyze_lmer rows: %d  runtime=%.2fs\n",
                nrow(lmer_res), lmer_runtime))
} else {
    cat("olink_lmer skipped\n")
}

# ---- Runtime record -------------------------------------------------------
writeLines(
    sprintf("OlinkAnalyze\tttest\t%.3fs\nOlinkAnalyze\tlmer\t%.3fs\n",
            runtime_s, lmer_runtime),
    file.path(out_dir, "olinkanalyze_runtime.tsv")
)
