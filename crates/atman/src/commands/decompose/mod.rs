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

/// Name of the search bound an automatic `k`-selection rule landed on,
/// or `null` when it resolved strictly inside its range.
///
/// A rule that returns its own ceiling has not selected anything — it
/// ran out of room, and the criterion it was asked to satisfy may never
/// have been met. `k_resolved: 30` with `k_max: 30` is indistinguishable
/// in a sidecar from a genuine selection that happens to land on 30, so
/// the distinction is recorded explicitly. Returns `None` for an
/// explicit `--k`, which is not a selection at all.
pub(super) fn k_selection_bound_hit(
    k_was_explicit: bool,
    k_resolved: usize,
    k_min: usize,
    k_max: usize,
) -> Option<&'static str> {
    if k_was_explicit || k_max <= k_min {
        return None;
    }
    if k_resolved >= k_max {
        Some("k_max")
    } else if k_resolved <= k_min {
        Some("k_min")
    } else {
        None
    }
}

/// Warn on stderr when an automatic `k`-selection rule resolved to one
/// of its own bounds. Silent otherwise.
pub(super) fn warn_if_k_selection_hit_bound(
    command: &str,
    rule: &str,
    k_resolved: usize,
    k_min: usize,
    k_max: usize,
) {
    match k_selection_bound_hit(false, k_resolved, k_min, k_max) {
        Some("k_max") => eprintln!(
            "{command}: warning: --k-selection {rule} resolved to k={k_resolved}, which is \
             --k-max. The criterion was not met anywhere in [{k_min}, {k_max}] — this is a \
             truncation, not a selection. Re-run with a larger --k-max to find the k the rule \
             actually wants."
        ),
        Some("k_min") => eprintln!(
            "{command}: warning: --k-selection {rule} resolved to k={k_resolved}, which is \
             --k-min. The rule was already satisfied at the smallest k searched; a smaller \
             --k-min may be appropriate."
        ),
        _ => {}
    }
}
