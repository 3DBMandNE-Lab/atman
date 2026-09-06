//! atman-core: platform-agnostic primitives and algorithms for proteomics data.
//!
//! The crate operates on the canonical Atman TSV schema. Upstream proteomics
//! ingest is handled by adapters (see `adapters/` at the repo root); this
//! crate provides the in-memory record types, statistical primitives, and
//! analysis algorithms applied after ingest.

pub mod align;
pub mod align_bootstrap;
pub mod align_project;
pub mod bench_decompose;
pub mod compositional;
pub mod contrast;
pub mod de;
pub mod decompose_unmix;
pub mod deqms;
pub mod ensemble;
pub mod expr;
pub mod fold_change;
pub mod gsea;
pub mod harmonize;
pub mod ica;
pub mod ica_mnar;
pub mod ica_null;
pub mod limma;
pub mod matrix;
pub mod modules_discover;
pub mod msqrob;
pub mod multivariate_t;
pub mod network;
pub mod network_differential;
pub mod nmf;
pub mod qc;
pub mod rng;
pub mod singscore;
pub mod stats;
pub mod studentized_range;
pub mod types;
pub mod variance_decomposition;

pub use align_bootstrap::{
    align_bootstrap, resample_rows, BootstrapParams, BootstrapRow, CohortMatrix, Decomposition,
};
pub use compositional::{apply_transform, Transform};
pub use contrast::{
    auc_mann_whitney, cohen_d_ci, cohen_d_pooled, kruskal_wallis, percentile, spearman_with_p,
    t_ci_mean, zscore, KruskalWallis,
};
pub use de::{bh_fdr, paired_t, PairedTResult, SkipReason};
pub use deqms::{deqms_shrink, tricube_moving_average, DeqmsShrinkage};
pub use ensemble::{
    aggregate_per_protein, assign_grade, combine_stouffer, EnsembleGrade, EnsembleInput,
    EnsembleRow, GradeThresholds,
};
pub use expr::{evaluate_column, split_top_level, Expr};
pub use fold_change::{
    compute_log2_fc, Comparison, FoldChangeInput, FoldChangeOutput, FoldChangePanel,
};
pub use gsea::{gsea, GseaConfig, GseaResult};
pub use ica_null::{archetype_null, generate_null_matrix, ArchetypeNullRow, NullMode, NullParams};
pub use matrix::{WidePanel, WidePanelRow};
pub use msqrob::{fit_msqrob, squeeze_variance, MsqrobFit, MsqrobOutcome};
pub use network::{
    adjacency, betweenness_centrality, eigenvector_centrality, influence_scores,
    pairwise_similarity, AdjacencyPolicy, InfluenceRow, SimilarityMatrix, SimilarityMetric,
};
pub use network_differential::{
    edge_pairwise_differential, edge_summary_differential, module_rewiring,
    signed_pairwise_correlations, CohortCorrelations, CohortData, EdgePairwiseRow, EdgeSummaryRow,
    ModuleRewiringRow, SignedMetric,
};
pub use rng::{derive_sub_seed, SplitMix64, Xoshiro256pp};
pub use singscore::{singscore, SingscoreRow};
pub use types::{
    Abundance, AssayId, Batch, DetectionLimit, MeasurementRecord, PeptideIdentity,
    PeptideMeasurementRecord, Platform, ProteinIdentity, QcFlag, Sample,
};
pub use variance_decomposition::{
    decompose_archetype_variance, FactorRow, FixedFactor, VarianceRow,
};
