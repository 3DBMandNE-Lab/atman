mod counterfactual;
mod ica;
mod nmf;
mod null;
mod unmix;
mod variance;

use counterfactual::run_counterfactual;
use ica::run_ica;
use nmf::nmf_run;
use null::run_null;
use unmix::run_unmix;
use variance::run_variance;

pub use counterfactual::CounterfactualArgs;
pub use ica::{IcaArgs, MissingnessModel, StabilityMetric, TransformArg};
pub use nmf::NmfArgs;
pub use null::{NullArgs, NullModeArg};
pub use unmix::UnmixArgs;
pub use variance::VarianceArgs;

use anyhow::Result;
use clap::{Args as ClapArgs, Subcommand};

#[derive(ClapArgs, Debug)]
pub struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Multi-seed FastICA decomposition with seed-stability reporting.
    Ica(IcaArgs),
    /// Permutation-null calibration of archetype stability.
    Null(NullArgs),
    /// Partition archetype activation variance into fixed-effect,
    /// random-intercept, and residual components per archetype.
    Variance(VarianceArgs),
    /// Geometric compartmental unmixing (VCA endmember extraction
    /// + FCLS/UCLS abundance estimation).
    Unmix(UnmixArgs),
    /// Reconstruct the abundance matrix with one program's activation
    /// set to a target value (default 0). Produces a delta matrix
    /// showing each program's per-protein, per-sample contribution.
    Counterfactual(CounterfactualArgs),
    /// Multiplicative-updates NMF (Frobenius loss, Lee & Seung).
    Nmf(NmfArgs),
}

pub fn run(args: Args) -> Result<()> {
    match args.command {
        Command::Ica(args) => run_ica(args),
        Command::Null(args) => run_null(args),
        Command::Variance(args) => run_variance(args),
        Command::Unmix(args) => run_unmix(args),
        Command::Counterfactual(args) => run_counterfactual(args),
        Command::Nmf(args) => nmf_run(args),
    }
}

/// Shared program-label formatter used across decomposition methods.
pub(super) fn program_name(index: usize) -> String {
    format!("program_{:02}", index + 1)
}
