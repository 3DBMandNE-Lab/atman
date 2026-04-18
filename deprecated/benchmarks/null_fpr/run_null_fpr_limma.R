#!/usr/bin/env Rscript
# Floor 3 limma null-FPR driver.
#
# For each of B sign-flip permutations per contrast, flip each subject's
# (A, B) labels with p=0.5, fit limma's lmFit + eBayes on the paired
# design, and record rejection counts at nominal BH q ∈ {0.01, 0.05, 0.10}.
#
# Invoked by run_null_fpr.py with matched seed; keeps permutation scheme
# consistent across Atman and limma up to the Python vs R RNG boundary
# (we do not share RNG state across processes, but B is large enough that
# cross-tool FPR comparison is unaffected).

suppressPackageStartupMessages({
    library(limma)
    library(dplyr)
    library(readr)
    library(tidyr)
})

# Minimal CLI parser: pairs of `--key value` args → named list.
parse_cli <- function(args) {
    out <- list()
    i <- 1
    while (i <= length(args)) {
        key <- sub("^--", "", args[i])
        val <- args[i + 1]
        out[[key]] <- val
        i <- i + 2
    }
    out
}
opts <- parse_cli(commandArgs(trailingOnly = TRUE))
if (is.null(opts$`input-dir`) || is.null(opts$output)) {
    stop("need --input-dir and --output")
}
opts$`n-perm` <- as.integer(opts$`n-perm`)
opts$seed    <- as.integer(opts$seed)
if (is.null(opts$contrasts)) opts$contrasts <- "PT1-PR1,PT2-PR2,PT2-PT1,PR2-PR1"

set.seed(opts$seed)
Q_CUTOFFS <- c(0.01, 0.05, 0.10)

# Load qc_measurements.tsv + samples.tsv, pivot to wide NPX matrix.
meas <- read_tsv(file.path(opts$`input-dir`, "qc_measurements.tsv"),
                 show_col_types = FALSE)
samp <- read_tsv(file.path(opts$`input-dir`, "samples.tsv"),
                 show_col_types = FALSE)

# Keep only measurements that passed QC (atman flag `dropped_by_qc == 0`).
if ("dropped_by_qc" %in% names(meas)) {
    meas <- meas %>% filter(dropped_by_qc == 0)
}

contrasts <- strsplit(opts$contrasts, ",")[[1]]

rows <- list()
for (cname in contrasts) {
    pair <- strsplit(cname, "-")[[1]]
    cond_a <- pair[1]; cond_b <- pair[2]

    # Samples relevant to this contrast.
    samp_c <- samp %>% filter(condition %in% c(cond_a, cond_b))
    meas_c <- meas %>% inner_join(samp_c, by = "sample_id")
    # Wide matrix: proteins × samples.
    wide <- meas_c %>%
        select(gene_symbol, sample_id, abundance) %>%
        group_by(gene_symbol, sample_id) %>% summarise(abundance = mean(abundance), .groups = "drop") %>%
        pivot_wider(names_from = sample_id, values_from = abundance)
    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$gene_symbol
    # Drop rows with any NA — required for paired lmFit with ~subject + condition.
    keep_rows <- rowSums(is.na(mat)) == 0
    mat <- mat[keep_rows, , drop = FALSE]

    meta <- samp_c[match(colnames(mat), samp_c$sample_id), ]
    stopifnot(!any(is.na(meta$condition)))

    n_tested <- nrow(mat)
    subjects <- unique(meta$subject_id)

    cat(sprintf("[limma] %s n_proteins=%d n_subjects=%d\n", cname, n_tested,
                length(subjects)), file = stderr())

    for (b in seq_len(opts$`n-perm`)) {
        # Sign-flip per subject with p=0.5.
        flip <- setNames(runif(length(subjects)) < 0.5, subjects)
        meta_b <- meta
        swap_idx <- flip[meta_b$subject_id]
        # Original labels a (= cond_a) and b (= cond_b). Swap for flipped subjects.
        meta_b$condition <- ifelse(
            swap_idx,
            ifelse(meta_b$condition == cond_a, cond_b, cond_a),
            meta_b$condition
        )

        # Skip perms where sign-flip ended up one-class or the balanced design broke.
        cond_factor <- factor(meta_b$condition, levels = c(cond_a, cond_b))
        subj_factor <- factor(meta_b$subject_id)
        if (length(unique(cond_factor)) < 2) next
        # Design may be rank-deficient for some permutations; guard.
        design <- tryCatch(
            model.matrix(~ subj_factor + cond_factor),
            error = function(e) NULL
        )
        if (is.null(design) || qr(design)$rank < ncol(design)) next

        fit <- tryCatch(
            { f <- lmFit(mat, design); eBayes(f) },
            error = function(e) NULL
        )
        if (is.null(fit)) next
        coef_name <- grep("^cond_factor", colnames(design), value = TRUE)
        if (length(coef_name) != 1) next
        tt <- topTable(fit, coef = coef_name, number = Inf,
                       sort.by = "none", adjust.method = "BH")

        for (qc in Q_CUTOFFS) {
            n_rej <- sum(tt$adj.P.Val < qc, na.rm = TRUE)
            rows[[length(rows) + 1]] <- data.frame(
                tool = "limma", mode = "eBayes", contrast = cname,
                perm = b - 1L, q_cutoff = qc, n_rejected = as.integer(n_rej),
                n_tested = n_tested
            )
        }
        if (b %% 25 == 0) {
            cat(sprintf("[limma] %s perm %d/%d\n", cname, b, opts$`n-perm`),
                file = stderr())
        }
    }
}

out <- bind_rows(rows)
write_tsv(out, opts$output)
cat(sprintf("[limma] wrote %s  (%d rows)\n", opts$output, nrow(out)),
    file = stderr())
