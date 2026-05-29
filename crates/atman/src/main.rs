//! atman: deterministic proteomics analysis engine.
//!
//! Atman operates on the canonical Atman TSV schema. Upstream proteomics
//! formats (Olink Explore NPX, SomaScan, MaxQuant/LFQ, DIA-NN, Spectronaut,
//! and any wide abundance matrix) land on the schema via adapters in
//! `adapters/`, or via the built-in `ingest-matrix` for stock wide-format
//! exports.

use anyhow::Result;
use clap::{Parser, Subcommand};

use atman::commands;

#[derive(Parser, Debug)]
#[command(
    name = "atman",
    version,
    about = "Deterministic proteomics analysis in Rust — Olink, SomaScan, DIA-NN, Spectronaut, MaxQuant"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Convert a wide proteomics matrix plus metadata into canonical Atman TSVs.
    IngestMatrix(Box<commands::ingest_matrix::Args>),
    /// Apply QC rule (`mask-warn-fail`: mark dropped_by_qc when qc_sample or qc_assay ≠ Pass).
    Qc(commands::qc::Args),
    /// Pivot QC'd long TSV into per-panel Dube-wide CSVs.
    Matrix(commands::matrix::Args),
    /// Cross-cohort program alignment with optional sensitivity sweep.
    Align(commands::align::Args),
    /// Multi-seed ICA decomposition with stability reporting.
    Decompose(commands::decompose::Args),
    /// Joint abundance-detection analysis for censored proteomics measurements.
    Detectability(commands::detectability::Args),
    /// Combine DE results across cohorts.
    Meta(commands::meta::Args),
    /// Compute log2 fold change per panel for the given comparisons.
    FoldChange(commands::fold_change::Args),
    /// Differential abundance with paired-t or moderated variance-shrinkage model.
    De(Box<commands::de::Args>),
    /// Gene-set enrichment analyses over DE results.
    Enrich(commands::enrich::Args),
    /// Compare matched contrast pairs for provocation-dependent asymmetry.
    Asymmetry(commands::asymmetry::Args),
    /// Bootstrap uncertainty for protein effects.
    Bootstrap(commands::bootstrap::Args),
    /// Subject-level program coupling with sign-consistency reporting.
    Coupling(commands::coupling::Args),
    /// Program-level utilities for ICA-derived program tables.
    Programs(commands::programs::Args),
    /// Compute LOO robustness and ranking stability from DE outputs.
    Robustness(commands::robustness::Args),
    /// Execute a pre-registered analysis plan (YAML/JSON) and emit a manifest.
    Run(commands::run::Args),
    /// Compute module trajectory scores from per-subject delta table and modules.tsv.
    ModuleTrajectory(commands::module_trajectory::Args),
    /// Module-level differential abundance (aggregate proteins into modules, test at module level).
    ModuleDe(commands::module_de::Args),
    /// Data-driven module discovery (WGCNA-style soft-threshold TOM + UPGMA).
    Modules(commands::modules::Args),
    /// Cross-tool benchmark harness for decomposition methods.
    Bench(commands::bench::Args),
    /// Null calibration by label permutation or paired sign flip.
    Null(commands::null::Args),
    /// Feature-covariance network analysis (`influence` subcommand).
    Network(commands::network::Args),
    /// Test subject-level log-ratios between two protein or module classes.
    Ratio(commands::ratio::Args),
    /// Recover TMT plex assignment from condition-invariant absence pattern.
    RecoverPlex(commands::recover_plex::Args),
    /// Cluster proteins by absence pattern (protein-side dual of recover-plex).
    AbsenceTopology(commands::absence_topology::Args),
    /// Leave-one-pair-out robustness diagnostic for paired DE.
    RobustPaired(commands::robust_paired::Args),
    /// Validate canonical Atman TSV inputs.
    Validate(commands::validate::Args),
    /// Generate compact reports from canonical Atman outputs.
    Report(commands::report::Args),
    /// Export nuisance-covariate-adjusted residual matrices.
    Residuals(commands::residuals::Args),
    /// Score biological modules from canonical Atman measurements.
    Score(commands::score::Args),
    /// Rank-transform measurements within each cohort/input directory.
    WithinCohortRank(commands::within_cohort_rank::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::IngestMatrix(args) => commands::ingest_matrix::run(*args),
        Command::Qc(args) => commands::qc::run(args),
        Command::Matrix(args) => commands::matrix::run(args),
        Command::Align(args) => commands::align::run(args),
        Command::Decompose(args) => commands::decompose::run(args),
        Command::Detectability(args) => commands::detectability::run(args),
        Command::Meta(args) => commands::meta::run(args),
        Command::FoldChange(args) => commands::fold_change::run(args),
        Command::De(args) => commands::de::run(*args),
        Command::Enrich(args) => commands::enrich::run(args),
        Command::Asymmetry(args) => commands::asymmetry::run(args),
        Command::Bootstrap(args) => commands::bootstrap::run(args),
        Command::Coupling(args) => commands::coupling::run(args),
        Command::Programs(args) => commands::programs::run(args),
        Command::Robustness(args) => commands::robustness::run(args),
        Command::Run(args) => commands::run::run(args),
        Command::ModuleTrajectory(args) => commands::module_trajectory::run(args),
        Command::ModuleDe(args) => commands::module_de::run(args),
        Command::Modules(args) => commands::modules::run(args),
        Command::Bench(args) => commands::bench::run(args),
        Command::Null(args) => commands::null::run(args),
        Command::Network(args) => commands::network::run(args),
        Command::Ratio(args) => commands::ratio::run(args),
        Command::RecoverPlex(args) => commands::recover_plex::run(args),
        Command::AbsenceTopology(args) => commands::absence_topology::run(args),
        Command::RobustPaired(args) => commands::robust_paired::run(args),
        Command::Validate(args) => commands::validate::run(args),
        Command::Report(args) => commands::report::run(args),
        Command::Residuals(args) => commands::residuals::run(args),
        Command::Score(args) => commands::score::run(args),
        Command::WithinCohortRank(args) => commands::within_cohort_rank::run(args),
    }
}
