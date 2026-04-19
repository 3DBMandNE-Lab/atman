//! atman: command-line engine for Olink Explore proteomics workflows.

use anyhow::Result;
use clap::{Parser, Subcommand};

use atman::commands;

#[derive(Parser, Debug)]
#[command(
    name = "atman",
    version,
    about = "Proteomics engine for Olink Explore NGS"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Parse raw NPX files into canonical long TSV + catalog + sample sheet.
    Ingest(commands::ingest::Args),
    /// Convert a wide proteomics matrix plus metadata into canonical Atman TSVs.
    IngestMatrix(Box<commands::ingest_matrix::Args>),
    /// Apply QC rule (Dube: mask rows where QC_Warning or Assay_Warning ≠ PASS).
    Qc(commands::qc::Args),
    /// Pivot QC'd long TSV into per-panel Dube-wide CSVs.
    Matrix(commands::matrix::Args),
    /// Combine DE results across cohorts.
    Meta(commands::meta::Args),
    /// Compute log2 fold change per panel for the given comparisons.
    FoldChange(commands::fold_change::Args),
    /// Differential abundance with paired-t or moderated variance-shrinkage model.
    De(commands::de::Args),
    /// Gene-set enrichment analyses over DE results.
    Enrich(commands::enrich::Args),
    /// Compare matched contrast pairs for provocation-dependent asymmetry.
    Asymmetry(commands::asymmetry::Args),
    /// Bootstrap uncertainty for protein effects.
    Bootstrap(commands::bootstrap::Args),
    /// Compute LOO robustness and ranking stability from DE outputs.
    Robustness(commands::robustness::Args),
    /// Compute module trajectory scores from per-subject delta table and modules.tsv.
    ModuleTrajectory(commands::module_trajectory::Args),
    /// Module-level differential abundance (aggregate proteins into modules, test at module level).
    ModuleDe(commands::module_de::Args),
    /// Null calibration by label permutation or paired sign flip.
    Null(commands::null::Args),
    /// Validate canonical Atman TSV inputs.
    Validate(commands::validate::Args),
    /// Generate compact reports from canonical Atman outputs.
    Report(commands::report::Args),
    /// Score biological modules from canonical Atman measurements.
    Score(commands::score::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Ingest(args) => commands::ingest::run(args),
        Command::IngestMatrix(args) => commands::ingest_matrix::run(*args),
        Command::Qc(args) => commands::qc::run(args),
        Command::Matrix(args) => commands::matrix::run(args),
        Command::Meta(args) => commands::meta::run(args),
        Command::FoldChange(args) => commands::fold_change::run(args),
        Command::De(args) => commands::de::run(args),
        Command::Enrich(args) => commands::enrich::run(args),
        Command::Asymmetry(args) => commands::asymmetry::run(args),
        Command::Bootstrap(args) => commands::bootstrap::run(args),
        Command::Robustness(args) => commands::robustness::run(args),
        Command::ModuleTrajectory(args) => commands::module_trajectory::run(args),
        Command::ModuleDe(args) => commands::module_de::run(args),
        Command::Null(args) => commands::null::run(args),
        Command::Validate(args) => commands::validate::run(args),
        Command::Report(args) => commands::report::run(args),
        Command::Score(args) => commands::score::run(args),
    }
}
