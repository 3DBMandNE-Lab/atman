#!/usr/bin/env Rscript
# limma on the Dube NPX matrix.
#
# Build wide NPX matrix (proteins x samples), fit a paired design (~ subject
# + condition) per contrast, write results in the shared schema.

suppressPackageStartupMessages({
    library(limma)
    library(dplyr)
    library(readr)
    library(tidyr)
})

if (!file.exists("Cargo.toml")) {
    stop("run from repo root; `Cargo.toml` not found in getwd() = ", getwd())
}

out_dir <- "benchmarks/out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load Dube NPX long table via OlinkAnalyze reader, then pivot ---------
suppressPackageStartupMessages(library(OlinkAnalyze))
npx_files <- c(
    "example_data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv",
    "example_data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv"
)
npx_raw <- bind_rows(lapply(npx_files, read_NPX))

parse_dube <- function(sid) {
    m <- regmatches(sid, regexec("^SSNA-([0-9]+)-(PR1|PT1|PR2|PT2)$", sid))
    do.call(rbind, lapply(m, function(x) if (length(x) == 3) x[2:3] else c(NA, NA)))
}
ids <- parse_dube(npx_raw$SampleID)
npx <- npx_raw %>%
    mutate(participant = ids[, 1], condition = ids[, 2]) %>%
    filter(!is.na(condition), Assay_Warning == "PASS", QC_Warning == "PASS")

# Wide: rows = OlinkID (unique per assay; a gene can appear on multiple
# panels — e.g. cardiometabolic + cardiometabolic_ii — so Assay alone is
# not a unique key). Gene symbol mapped back at write time.
olink_to_gene <- npx %>% select(OlinkID, Assay) %>% distinct()
wide <- npx %>% select(OlinkID, SampleID, NPX) %>%
    pivot_wider(names_from = SampleID, values_from = NPX)
mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$OlinkID
# Keep rows with at least one non-NA value — limma handles per-row NAs.
mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

# Rebuild sample metadata in column order.
meta <- npx %>% select(SampleID, participant, condition) %>% distinct()
meta <- meta[match(colnames(mat), meta$SampleID), ]
stopifnot(!any(is.na(meta$condition)))

contrasts_def <- list(
    "PT1-PR1" = c("PT1", "PR1"),
    "PT2-PR2" = c("PT2", "PR2"),
    "PT2-PT1" = c("PT2", "PT1"),
    "PR2-PR1" = c("PR2", "PR1")
)

all_de <- list()
t_start <- Sys.time()
for (cname in names(contrasts_def)) {
    pair <- contrasts_def[[cname]]
    keep <- meta$condition %in% pair
    sub_mat <- mat[, keep, drop = FALSE]
    sub_meta <- meta[keep, ]
    sub_meta$condition <- factor(sub_meta$condition, levels = pair)
    sub_meta$participant <- factor(sub_meta$participant)

    # Paired design: block = participant, coef = condition.
    design <- model.matrix(~ participant + condition, data = sub_meta)
    fit <- lmFit(sub_mat, design)
    fit <- eBayes(fit)
    coef_name <- grep("^condition", colnames(design), value = TRUE)
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none",
                   adjust.method = "BH")
    tt$OlinkID <- rownames(tt)
    tt <- left_join(tt, olink_to_gene, by = "OlinkID")
    tt$gene <- tt$Assay
    tt$contrast <- cname
    all_de[[cname]] <- tt
}
t_end <- Sys.time()
runtime_s <- as.numeric(difftime(t_end, t_start, units = "secs"))

de <- bind_rows(all_de) %>%
    transmute(
        contrast  = contrast,
        gene      = gene,
        mean_diff = logFC,
        pvalue    = P.Value,
        qvalue    = adj.P.Val,
        tool      = "limma",
        runtime_s = runtime_s,
        peak_rss_mb = NA_real_
    )
write_tsv(de, file.path(out_dir, "limma_de.tsv"))
cat(sprintf("limma DE rows: %d  runtime=%.2fs\n", nrow(de), runtime_s))
