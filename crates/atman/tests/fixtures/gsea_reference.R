# Reference enrichment scores from fgsea::fgseaSimple on a planted
# fixture. Atman's `enrich gsea` must reproduce the ES values within
# tight numerical tolerance. P-values and NES depend on the
# permutation null draws and differ across PRNG implementations;
# the parity claim is on the deterministic ES.
#
# Run: Rscript gsea_reference.R
# (Set ATMAN_FIXTURE_DIR to override the output directory.)

suppressPackageStartupMessages({
  library(fgsea)
})

script_arg <- Sys.getenv("ATMAN_FIXTURE_DIR")
fixture_dir <- if (nzchar(script_arg)) script_arg else {
  args <- commandArgs(trailingOnly = FALSE)
  self <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(self) == 1 && file.exists(self)) dirname(normalizePath(self))
  else "crates/atman/tests/fixtures"
}

set.seed(20260428L)

# 300 genes ranked descending by stat. The top 12 ("HIT001"..."HIT012")
# are clear positives; the bottom 8 ("DOWN01"..."DOWN08") are clear
# negatives; the remaining 280 are noise centered near zero. Three
# planted gene sets:
#   * "top_loaded" — exactly the 12 HIT genes (ES → +1)
#   * "bottom_loaded" — exactly the 8 DOWN genes (ES → -1)
#   * "scattered" — every 30th noise gene (ES near 0)

n_hit <- 12L
n_down <- 8L
n_noise <- 280L
hit_names <- sprintf("HIT%03d", seq_len(n_hit))
down_names <- sprintf("DOWN%02d", seq_len(n_down))
noise_names <- sprintf("NOISE%03d", seq_len(n_noise))

# Construct stats by descending order so the ranked list is also the
# input gene order (no re-sorting needed downstream).
hit_stats <- seq(from = 5.0, to = 4.0, length.out = n_hit)
noise_stats <- sort(rnorm(n_noise, mean = 0, sd = 0.5), decreasing = TRUE)
down_stats <- seq(from = -4.0, to = -5.0, length.out = n_down)
ranked_genes <- c(hit_names, noise_names, down_names)
ranked_stats <- c(hit_stats, noise_stats, down_stats)
stopifnot(!is.unsorted(rev(ranked_stats)))  # descending

# Emit the ranked list as a TSV in atman de_results.tsv shape so the
# integration test can feed it directly to `atman enrich gsea`.
de <- data.frame(
  panel = "p",
  assay_id = ranked_genes,
  gene_symbol = ranked_genes,
  uniprot = "",
  comparison = "Case-Control",
  n_pairs = 10L,
  mean_a = "",
  mean_b = "",
  mean_diff = ranked_stats,
  t = ranked_stats,
  df = 8L,
  # `signed_log10_p` ranking should reproduce ranked_stats. To keep
  # this faithful for atman's default rank-by, encode p as
  # 10^(-|stat|) so sign(mean_diff) * -log10(p) == stat.
  p_value = 10^(-abs(ranked_stats)),
  bh_q = 10^(-abs(ranked_stats)),
  skip_reason = "",
  stringsAsFactors = FALSE
)
write.table(
  de,
  file = file.path(fixture_dir, "gsea_de_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# Gene sets in atman gene_sets.tsv shape (long form: set_name, gene_symbol).
sets_long <- rbind(
  data.frame(set_name = "top_loaded", gene_symbol = hit_names, stringsAsFactors = FALSE),
  data.frame(set_name = "bottom_loaded", gene_symbol = down_names, stringsAsFactors = FALSE),
  data.frame(
    set_name = "scattered",
    gene_symbol = noise_names[seq(1L, n_noise, by = 30L)],
    stringsAsFactors = FALSE
  )
)
write.table(
  sets_long,
  file = file.path(fixture_dir, "gsea_gene_sets.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# Compute the reference ES with fgseaSimple. minSize=1 so the small
# scattered set isn't filtered.
fgsea_input <- setNames(ranked_stats, ranked_genes)
fgsea_pathways <- list(
  top_loaded = hit_names,
  bottom_loaded = down_names,
  scattered = noise_names[seq(1L, n_noise, by = 30L)]
)
res <- fgseaSimple(
  pathways = fgsea_pathways,
  stats = fgsea_input,
  nperm = 1000L,
  minSize = 1L,
  maxSize = 500L
)
ref <- data.frame(
  set_name = as.character(res$pathway),
  set_size = as.integer(res$size),
  es = as.numeric(res$ES),
  stringsAsFactors = FALSE
)
ref <- ref[order(ref$set_name), ]
out_path <- file.path(fixture_dir, "gsea_reference.tsv")
write.table(ref, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", out_path, "with", nrow(ref), "rows\n")
print(ref)
