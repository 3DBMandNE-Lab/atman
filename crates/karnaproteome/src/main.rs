//! karnaproteome: engine binary for proteomics data. v0.1 ships four subcommands
//! targeting reproduction of the Dube et al. Olink Explore published outputs.

use anyhow::Result;
use clap::{Parser, Subcommand};

use karnaproteome::commands;

#[derive(Parser, Debug)]
#[command(
    name = "karnaproteome",
    version,
    about = "Proteomics engine (Olink Explore NGS, v0.1)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Parse raw NPX files into canonical long TSV + catalog + sample sheet.
    Ingest(commands::ingest::Args),
    /// Apply QC rule (Dube: mask rows where QC_Warning or Assay_Warning ≠ PASS).
    Qc(commands::qc::Args),
    /// Pivot QC'd long TSV into per-panel Dube-wide CSVs.
    Matrix(commands::matrix::Args),
    /// Compute log2 fold change per panel for the given comparisons.
    FoldChange(commands::fold_change::Args),
    /// Differential abundance with paired-t or moderated variance-shrinkage model.
    De(commands::de::Args),
    /// Compare matched contrast pairs for provocation-dependent asymmetry.
    Asymmetry(commands::asymmetry::Args),
    /// Compute LOO robustness and ranking stability from DE outputs.
    Robustness(commands::robustness::Args),
    /// Compute module trajectory scores from per-subject delta table and modules.tsv.
    ModuleTrajectory(commands::module_trajectory::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Ingest(args) => commands::ingest::run(args),
        Command::Qc(args) => commands::qc::run(args),
        Command::Matrix(args) => commands::matrix::run(args),
        Command::FoldChange(args) => commands::fold_change::run(args),
        Command::De(args) => commands::de::run(args),
        Command::Asymmetry(args) => commands::asymmetry::run(args),
        Command::Robustness(args) => commands::robustness::run(args),
        Command::ModuleTrajectory(args) => commands::module_trajectory::run(args),
    }
}
