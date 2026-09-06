//! Cross-cohort harmonisation with a fit-then-apply contract.
//!
//! The point of the contract is that it makes leakage structurally
//! impossible rather than a discipline: [`fit`] sees the training
//! cohorts and writes a [`HarmonizeModel`]; [`apply`] sees that model
//! and one new cohort, and nothing else. A method that cannot express
//! itself as a transportable model cannot be benchmarked this way, and
//! that is a fact about the method rather than a limitation of the
//! harness — location-scale correctors estimate every batch's
//! parameters jointly, so there is no artefact to carry to a cohort
//! that was not in the fit.
//!
//! The model carries two things: how to map an incoming cohort onto the
//! shared feature space, and a disease direction learned on the training
//! cohorts in that method's own representation. Applying it to a
//! held-out cohort yields one transfer score per subject, which is what
//! gets tested there.

use std::collections::BTreeMap;

use crate::rng::SplitMix64;

/// How a method maps a cohort onto the shared representation.
#[derive(Debug, Clone, PartialEq)]
pub enum HarmonizeMethod {
    /// Per-cohort standardisation. The null harmonisation: it removes
    /// cohort location and scale and nothing else. Requires no fitted
    /// state, since a new cohort is standardised by its own statistics.
    ZScore,
    /// Within-sample rank, rescaled to `[0, 1]`. Also stateless, so it
    /// is trivially appliable — worth knowing when a comparator scores
    /// well that it learned nothing to do so.
    Rank,
    /// Quantile normalisation against a reference profile learned on the
    /// training cohorts. Unlike rank this has real fitted state.
    Quantile { reference: Vec<f64> },
    /// Per-sample division by the geometric mean of a reference protein
    /// set, in the style of reference-protein normalisation for CSF.
    ///
    /// The reference set is chosen on the training cohorts only, which
    /// the fit/apply split enforces: [`apply`] receives the chosen
    /// indices in the model and has no way to reselect.
    ReferenceProtein { reference_features: Vec<String> },
}

impl HarmonizeMethod {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::ZScore => "zscore",
            Self::Rank => "rank",
            Self::Quantile { .. } => "quantile",
            Self::ReferenceProtein { .. } => "reference-protein",
        }
    }

    /// Whether the method learns anything from the training cohorts.
    ///
    /// A stateless method scoring as well as a fitted one is worth
    /// reporting rather than hiding: it means the fitted state was not
    /// what carried the signal.
    pub fn is_fitted(&self) -> bool {
        !matches!(self, Self::ZScore | Self::Rank)
    }
}

/// One cohort's subject × feature matrix on a named feature axis.
#[derive(Debug, Clone)]
pub struct CohortData {
    pub label: String,
    pub subject_ids: Vec<String>,
    pub features: Vec<String>,
    /// `[subject][feature]`; non-finite means unobserved.
    pub values: Vec<Vec<f64>>,
    /// `true` for a case, `false` for a control.
    pub is_case: Vec<bool>,
}

/// A transportable harmonisation model plus the direction learned under
/// it. Written by [`fit`], consumed by [`apply`].
#[derive(Debug, Clone)]
pub struct HarmonizeModel {
    pub method: HarmonizeMethod,
    /// Shared feature axis: the intersection over training cohorts.
    pub features: Vec<String>,
    /// Per-feature disease direction in the harmonised representation.
    pub direction: Vec<f64>,
    /// Labels of the cohorts this was fit on, so a replay can prove a
    /// held-out cohort was absent.
    pub fit_cohorts: Vec<String>,
    /// Whether the direction was learned on shuffled labels.
    ///
    /// The permuted arm is the negative control: every method scores
    /// something on a held-out cohort, so a benchmark where every entry
    /// beats zero measures nothing. A method whose permuted arm still
    /// separates cases is detecting structure unrelated to disease.
    pub permuted: bool,
}

fn mean_sd(xs: &[f64]) -> (f64, f64) {
    let obs: Vec<f64> = xs.iter().copied().filter(|v| v.is_finite()).collect();
    if obs.len() < 2 {
        return (0.0, 0.0);
    }
    let mean = obs.iter().sum::<f64>() / obs.len() as f64;
    let var = obs.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (obs.len() as f64 - 1.0);
    (mean, var.sqrt())
}

/// Within-sample ranks scaled to `[0, 1]`, ties averaged. Unobserved
/// features stay non-finite rather than being ranked as if present.
fn rank_row(row: &[f64]) -> Vec<f64> {
    let obs: Vec<(usize, f64)> = row
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite())
        .map(|(i, v)| (i, *v))
        .collect();
    let mut out = vec![f64::NAN; row.len()];
    if obs.is_empty() {
        return out;
    }
    let mut order = obs.clone();
    order.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let n = order.len() as f64;
    let mut i = 0usize;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && order[j + 1].1 == order[i].1 {
            j += 1;
        }
        let avg_rank = ((i + j) as f64) / 2.0;
        for item in order.iter().take(j + 1).skip(i) {
            out[item.0] = if n > 1.0 { avg_rank / (n - 1.0) } else { 0.5 };
        }
        i = j + 1;
    }
    out
}

/// Map one cohort into the model's representation on the model's
/// feature axis. Features the cohort lacks stay non-finite.
fn represent(
    model_features: &[String],
    method: &HarmonizeMethod,
    cohort: &CohortData,
) -> Vec<Vec<f64>> {
    let index: BTreeMap<&str, usize> = cohort
        .features
        .iter()
        .enumerate()
        .map(|(i, f)| (f.as_str(), i))
        .collect();
    // Restrict to the model's axis first, so every method sees the same
    // features in the same order.
    let restricted: Vec<Vec<f64>> = cohort
        .values
        .iter()
        .map(|row| {
            model_features
                .iter()
                .map(|f| index.get(f.as_str()).map(|j| row[*j]).unwrap_or(f64::NAN))
                .collect()
        })
        .collect();

    match method {
        HarmonizeMethod::ZScore => {
            let p = model_features.len();
            let stats: Vec<(f64, f64)> = (0..p)
                .map(|j| {
                    let col: Vec<f64> = restricted.iter().map(|r| r[j]).collect();
                    mean_sd(&col)
                })
                .collect();
            restricted
                .iter()
                .map(|row| {
                    row.iter()
                        .zip(stats.iter())
                        .map(|(v, (m, s))| if *s > 0.0 { (v - m) / s } else { f64::NAN })
                        .collect()
                })
                .collect()
        }
        HarmonizeMethod::Rank => restricted.iter().map(|r| rank_row(r)).collect(),
        HarmonizeMethod::Quantile { reference } => restricted
            .iter()
            .map(|row| {
                let ranks = rank_row(row);
                ranks
                    .iter()
                    .map(|r| {
                        if !r.is_finite() || reference.is_empty() {
                            return f64::NAN;
                        }
                        // Map the within-sample quantile onto the
                        // training-cohort reference profile.
                        let pos = r * (reference.len() as f64 - 1.0);
                        let lo = pos.floor() as usize;
                        let hi = pos.ceil() as usize;
                        let frac = pos - lo as f64;
                        reference[lo] * (1.0 - frac) + reference[hi.min(reference.len() - 1)] * frac
                    })
                    .collect()
            })
            .collect(),
        HarmonizeMethod::ReferenceProtein { reference_features } => {
            let ref_idx: Vec<usize> = reference_features
                .iter()
                .filter_map(|f| model_features.iter().position(|g| g == f))
                .collect();
            restricted
                .iter()
                .map(|row| {
                    // Geometric mean on a log scale is the arithmetic
                    // mean; the inputs here are already log abundances.
                    let refs: Vec<f64> = ref_idx
                        .iter()
                        .map(|j| row[*j])
                        .filter(|v| v.is_finite())
                        .collect();
                    if refs.is_empty() {
                        return vec![f64::NAN; row.len()];
                    }
                    let anchor = refs.iter().sum::<f64>() / refs.len() as f64;
                    row.iter().map(|v| v - anchor).collect()
                })
                .collect()
        }
    }
}

/// Feature intersection across cohorts, in stable order.
fn shared_features(cohorts: &[CohortData]) -> Vec<String> {
    let Some(first) = cohorts.first() else {
        return Vec::new();
    };
    first
        .features
        .iter()
        .filter(|f| {
            cohorts
                .iter()
                .skip(1)
                .all(|c| c.features.iter().any(|g| g == *f))
        })
        .cloned()
        .collect()
}

/// Fit a harmonisation model and a disease direction on the training
/// cohorts.
///
/// `permute` shuffles case/control labels within each cohort before
/// learning the direction, preserving group sizes. That is the negative
/// control arm and it travels through the identical apply path.
pub fn fit(
    cohorts: &[CohortData],
    method_spec: MethodSpec,
    permute: bool,
    seed: u64,
) -> Result<HarmonizeModel, String> {
    if cohorts.len() < 2 {
        return Err("harmonize fit: at least 2 training cohorts required".into());
    }
    let features = shared_features(cohorts);
    if features.is_empty() {
        return Err("harmonize fit: training cohorts share no features".into());
    }

    // Resolve any state the method learns from the training data.
    let method = match method_spec {
        MethodSpec::ZScore => HarmonizeMethod::ZScore,
        MethodSpec::Rank => HarmonizeMethod::Rank,
        MethodSpec::Quantile => {
            // Reference profile: the mean sorted observed vector across
            // all training subjects, on the shared axis.
            let mut sorted_rows: Vec<Vec<f64>> = Vec::new();
            for c in cohorts {
                let rep = represent(&features, &HarmonizeMethod::ZScore, c);
                for row in rep {
                    let mut obs: Vec<f64> = row.into_iter().filter(|v| v.is_finite()).collect();
                    if obs.len() == features.len() {
                        obs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                        sorted_rows.push(obs);
                    }
                }
            }
            if sorted_rows.is_empty() {
                return Err(
                    "harmonize fit: quantile needs at least one training subject observed on \
                     every shared feature"
                        .into(),
                );
            }
            let reference: Vec<f64> = (0..features.len())
                .map(|j| sorted_rows.iter().map(|r| r[j]).sum::<f64>() / sorted_rows.len() as f64)
                .collect();
            HarmonizeMethod::Quantile { reference }
        }
        MethodSpec::ReferenceProtein { k } => {
            // Smallest within-cohort variance of log abundance, averaged
            // across training cohorts, restricted to features observed
            // in every training cohort.
            let mut scores: Vec<(f64, &String)> = Vec::new();
            for (j, f) in features.iter().enumerate() {
                let mut per_cohort: Vec<f64> = Vec::new();
                let mut observed_everywhere = true;
                for c in cohorts {
                    let idx = c.features.iter().position(|g| g == f);
                    let Some(idx) = idx else {
                        observed_everywhere = false;
                        break;
                    };
                    let col: Vec<f64> = c.values.iter().map(|r| r[idx]).collect();
                    let (_, sd) = mean_sd(&col);
                    if sd <= 0.0 {
                        observed_everywhere = false;
                        break;
                    }
                    per_cohort.push(sd * sd);
                }
                let _ = j;
                if observed_everywhere && !per_cohort.is_empty() {
                    scores.push((per_cohort.iter().sum::<f64>() / per_cohort.len() as f64, f));
                }
            }
            if scores.len() < k {
                return Err(format!(
                    "harmonize fit: only {} features are observed with non-zero variance in \
                     every training cohort; --reference-k is {k}",
                    scores.len()
                ));
            }
            scores.sort_by(|a, b| {
                a.0.partial_cmp(&b.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.1.cmp(b.1))
            });
            HarmonizeMethod::ReferenceProtein {
                reference_features: scores.iter().take(k).map(|(_, f)| (*f).clone()).collect(),
            }
        }
    };

    // Learn the direction in this method's representation.
    let mut case_sums = vec![0.0_f64; features.len()];
    let mut case_ns = vec![0usize; features.len()];
    let mut ctrl_sums = vec![0.0_f64; features.len()];
    let mut ctrl_ns = vec![0usize; features.len()];
    let mut rng = SplitMix64::new(seed);
    for c in cohorts {
        let rep = represent(&features, &method, c);
        let mut labels = c.is_case.clone();
        if permute {
            // Shuffle within cohort, preserving group sizes.
            for i in (1..labels.len()).rev() {
                let j = (rng.next_u64() % (i as u64 + 1)) as usize;
                labels.swap(i, j);
            }
        }
        for (row, is_case) in rep.iter().zip(labels.iter()) {
            for (j, v) in row.iter().enumerate() {
                if !v.is_finite() {
                    continue;
                }
                if *is_case {
                    case_sums[j] += v;
                    case_ns[j] += 1;
                } else {
                    ctrl_sums[j] += v;
                    ctrl_ns[j] += 1;
                }
            }
        }
    }
    let direction: Vec<f64> = (0..features.len())
        .map(|j| {
            if case_ns[j] == 0 || ctrl_ns[j] == 0 {
                0.0
            } else {
                case_sums[j] / case_ns[j] as f64 - ctrl_sums[j] / ctrl_ns[j] as f64
            }
        })
        .collect();

    Ok(HarmonizeModel {
        method,
        features,
        direction,
        fit_cohorts: cohorts.iter().map(|c| c.label.clone()).collect(),
        permuted: permute,
    })
}

/// Which method to fit, before its training-dependent state exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MethodSpec {
    ZScore,
    Rank,
    Quantile,
    ReferenceProtein { k: usize },
}

impl MethodSpec {
    pub fn parse(s: &str, reference_k: usize) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "zscore" => Some(Self::ZScore),
            "rank" => Some(Self::Rank),
            "quantile" => Some(Self::Quantile),
            "reference-protein" => Some(Self::ReferenceProtein { k: reference_k }),
            _ => None,
        }
    }
}

/// One held-out subject's transfer score.
#[derive(Debug, Clone)]
pub struct TransferScore {
    pub subject_id: String,
    pub score: f64,
    pub is_case: bool,
    /// Shared features this subject actually contributed.
    pub n_features_used: usize,
}

/// Apply a fitted model to a cohort it has never seen.
///
/// Takes the model and one cohort. It cannot see the training data, so
/// there is nothing to leak.
/// How the transfer score treats a direction's mean component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectionCentering {
    /// Score with the direction as learned.
    None,
    /// Subtract the direction's mean over each subject's own used
    /// features before scoring.
    ///
    /// Every per-sample normalisation leaves a per-sample term in the
    /// harmonised value, and a direction with a non-zero mean projects
    /// onto it. For reference-protein normalisation the algebra is
    /// exact: `x'_ij = x_ij − a_i`, so the score carries
    /// `−a_i · mean(w)`. If the anchor `a_i` differs between arms in the
    /// held-out cohort, ANY direction separates them, including one
    /// learned from shuffled labels. Centring the direction over each
    /// subject's used set makes that term identically zero.
    ///
    /// The cost is honest rather than hidden: a genuine disease effect
    /// that raises every protein uniformly is indistinguishable from a
    /// per-sample anchor shift under a per-sample normalisation, so
    /// centring does not discard recoverable signal, it declines to
    /// claim unrecoverable signal.
    PerSubject,
}

pub fn apply(model: &HarmonizeModel, cohort: &CohortData) -> Result<Vec<TransferScore>, String> {
    apply_with(model, cohort, DirectionCentering::None)
}

pub fn apply_with(
    model: &HarmonizeModel,
    cohort: &CohortData,
    centering: DirectionCentering,
) -> Result<Vec<TransferScore>, String> {
    if cohort.subject_ids.len() != cohort.values.len() {
        return Err("harmonize apply: subject count does not match matrix rows".into());
    }
    if model.fit_cohorts.contains(&cohort.label) {
        return Err(format!(
            "harmonize apply: cohort {:?} was one of the model's training cohorts ({:?}); \
             applying a model to data it was fit on is not a held-out evaluation",
            cohort.label, model.fit_cohorts
        ));
    }
    let rep = represent(&model.features, &model.method, cohort);
    Ok(rep
        .iter()
        .enumerate()
        .map(|(i, row)| {
            // Direction mean over THIS subject's used features. Which
            // features are usable varies by subject under MNAR
            // dropout, so a globally centred direction would still
            // leave a residual mean per subject.
            let mut used = 0usize;
            let mut w_sum = 0.0_f64;
            for (j, v) in row.iter().enumerate() {
                if v.is_finite() && model.direction[j] != 0.0 {
                    w_sum += model.direction[j];
                    used += 1;
                }
            }
            let shift = match centering {
                DirectionCentering::None => 0.0,
                DirectionCentering::PerSubject if used > 0 => w_sum / used as f64,
                DirectionCentering::PerSubject => 0.0,
            };
            let mut dot = 0.0_f64;
            for (j, v) in row.iter().enumerate() {
                if v.is_finite() && model.direction[j] != 0.0 {
                    dot += v * (model.direction[j] - shift);
                }
            }
            TransferScore {
                subject_id: cohort.subject_ids[i].clone(),
                score: if used == 0 {
                    f64::NAN
                } else {
                    dot / used as f64
                },
                is_case: cohort.is_case.get(i).copied().unwrap_or(false),
                n_features_used: used,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Planted-signal cohort. `Spec` keeps the argument count down and
    /// makes the signal fraction explicit, which turns out to matter:
    /// see [`the_permuted_arm_does_not_separate_the_held_out_cohort`].
    #[derive(Clone, Copy)]
    pub(crate) struct Spec {
        n: usize,
        p: usize,
        n_signal: usize,
        offset: f64,
        scale: f64,
        effect: f64,
    }

    pub(crate) fn cohort(label: &str, spec: Spec, seed: u64) -> CohortData {
        let mut rng = SplitMix64::new(seed);
        let mut values = Vec::with_capacity(spec.n);
        let mut is_case = Vec::with_capacity(spec.n);
        for i in 0..spec.n {
            let case = i % 2 == 0;
            is_case.push(case);
            let row: Vec<f64> = (0..spec.p)
                .map(|j| {
                    let noise = (rng.next_u64() % 1000) as f64 / 1000.0 - 0.5;
                    let signal = if j < spec.n_signal && case {
                        spec.effect
                    } else {
                        0.0
                    };
                    spec.offset + spec.scale * (noise + signal)
                })
                .collect();
            values.push(row);
        }
        CohortData {
            label: label.to_string(),
            subject_ids: (0..spec.n).map(|i| format!("{label}_S{i:03}")).collect(),
            features: (0..spec.p).map(|j| format!("F{j:03}")).collect(),
            values,
            is_case,
        }
    }

    /// 10 signal features in 200 is a realistic sparsity. An earlier
    /// version of this fixture used 12 in 60, and at that density a
    /// RANDOM direction aligns with the signal subspace often enough
    /// that the permuted null ran from −10.7 to +12.0 against a real
    /// effect of +11.9 — the negative control could not separate itself
    /// from the result. The arm was working; the fixture was not.
    pub(crate) const SPEC: Spec = Spec {
        n: 60,
        p: 200,
        n_signal: 10,
        offset: 10.0,
        scale: 1.0,
        effect: 1.5,
    };

    pub(crate) fn training() -> Vec<CohortData> {
        vec![
            cohort("a", SPEC, 1),
            cohort(
                "b",
                Spec {
                    offset: 25.0,
                    scale: 3.0,
                    ..SPEC
                },
                2,
            ),
        ]
    }

    pub(crate) fn held_out() -> CohortData {
        cohort(
            "c",
            Spec {
                offset: 100.0,
                scale: 7.0,
                ..SPEC
            },
            3,
        )
    }

    pub(crate) fn separation(scores: &[TransferScore]) -> f64 {
        let case: Vec<f64> = scores
            .iter()
            .filter(|s| s.is_case && s.score.is_finite())
            .map(|s| s.score)
            .collect();
        let ctrl: Vec<f64> = scores
            .iter()
            .filter(|s| !s.is_case && s.score.is_finite())
            .map(|s| s.score)
            .collect();
        let (mc, sc) = mean_sd(&case);
        let (mk, sk) = mean_sd(&ctrl);
        let pooled = ((sc * sc + sk * sk) / 2.0).sqrt();
        if pooled > 0.0 {
            (mc - mk) / pooled
        } else {
            0.0
        }
    }

    /// The contract: applying a model to a cohort it was fit on is not a
    /// held-out evaluation and must be refused by the interface rather
    /// than avoided by discipline.
    #[test]
    fn apply_refuses_a_training_cohort() {
        let train = training();
        let model = fit(&train, MethodSpec::ZScore, false, 42).unwrap();
        let err = apply(&model, &train[0]).unwrap_err();
        assert!(err.contains("training cohorts"), "unexpected error: {err}");
    }

    /// A real disease direction should transfer to a cohort with a
    /// different offset and scale.
    #[test]
    fn a_learned_direction_transfers_to_a_held_out_cohort() {
        let train = training();
        let held_out = held_out();
        for spec in [
            MethodSpec::ZScore,
            MethodSpec::Rank,
            MethodSpec::Quantile,
            MethodSpec::ReferenceProtein { k: 8 },
        ] {
            let model = fit(&train, spec, false, 42).unwrap();
            let scores = apply(&model, &held_out).unwrap();
            assert_eq!(scores.len(), SPEC.n);
            assert!(
                separation(&scores) > 0.5,
                "{:?} failed to transfer: separation {}",
                spec,
                separation(&scores)
            );
        }
    }

    /// The negative control. With labels shuffled in training there is
    /// no disease direction to carry, so the held-out separation must
    /// collapse. A method that still separates is detecting something
    /// other than disease, which is the whole point of the arm.
    /// The negative control, tested against a DISTRIBUTION rather than
    /// one draw.
    ///
    /// A single permuted fit is a single sample from the null, and an
    /// earlier version of this test asserted against exactly one. It
    /// passed, and it was luck: on the fixture then in use the null
    /// spanned the real effect, so the arm and the result were
    /// indistinguishable and the test said nothing. Every permuted draw
    /// must now fall below the real effect.
    #[test]
    fn the_permuted_arm_does_not_separate_the_held_out_cohort() {
        let train = training();
        let held_out = held_out();
        for spec in [MethodSpec::ZScore, MethodSpec::Rank, MethodSpec::Quantile] {
            let real_sep =
                separation(&apply(&fit(&train, spec, false, 42).unwrap(), &held_out).unwrap())
                    .abs();
            let null: Vec<f64> = (1..=12)
                .map(|seed| {
                    let m = fit(&train, spec, true, seed).unwrap();
                    assert!(m.permuted);
                    separation(&apply(&m, &held_out).unwrap()).abs()
                })
                .collect();
            let worst = null.iter().cloned().fold(0.0_f64, f64::max);
            assert!(
                worst < real_sep,
                "{spec:?}: the largest of 12 permuted draws ({worst:.2}) reached the real effect \
                 ({real_sep:.2}), so the control cannot distinguish itself from the result.\n\
                 null: {null:?}"
            );
        }
    }

    /// Reference proteins are selected on training cohorts only, and the
    /// chosen set travels in the model. `apply` has no selection code.
    #[test]
    fn reference_proteins_are_chosen_on_training_data_and_carried() {
        let train = training();
        let model = fit(&train, MethodSpec::ReferenceProtein { k: 8 }, false, 42).unwrap();
        match &model.method {
            HarmonizeMethod::ReferenceProtein { reference_features } => {
                assert_eq!(reference_features.len(), 8);
                // Low-variance features are the non-signal ones here.
                assert!(
                    reference_features.iter().all(|f| {
                        let idx: usize = f.trim_start_matches('F').parse().unwrap();
                        idx >= 12
                    }),
                    "selection should avoid the planted signal features: {reference_features:?}"
                );
            }
            other => panic!("wrong method resolved: {other:?}"),
        }
    }

    /// Stateless methods must be labelled as such: a stateless method
    /// scoring well means the fitted state was not what carried it.
    #[test]
    fn stateless_methods_are_identified() {
        let train = training();
        assert!(!fit(&train, MethodSpec::ZScore, false, 1)
            .unwrap()
            .method
            .is_fitted());
        assert!(!fit(&train, MethodSpec::Rank, false, 1)
            .unwrap()
            .method
            .is_fitted());
        assert!(fit(&train, MethodSpec::Quantile, false, 1)
            .unwrap()
            .method
            .is_fitted());
        assert!(fit(&train, MethodSpec::ReferenceProtein { k: 4 }, false, 1)
            .unwrap()
            .method
            .is_fitted());
    }

    #[test]
    fn fit_refuses_a_single_cohort_and_disjoint_features() {
        let one = vec![cohort("a", SPEC, 1)];
        assert!(fit(&one, MethodSpec::ZScore, false, 1).is_err());

        let mut b = cohort("b", SPEC, 2);
        b.features = (0..SPEC.p).map(|j| format!("OTHER{j:03}")).collect();
        let disjoint = vec![cohort("a", SPEC, 1), b];
        let err = fit(&disjoint, MethodSpec::ZScore, false, 1).unwrap_err();
        assert!(err.contains("share no features"), "unexpected: {err}");
    }

    /// Reproduces the artefact the CSF session measured: under
    /// reference-protein normalisation, a held-out cohort whose anchor
    /// differs between arms is separated by a direction learned from
    /// SHUFFLED labels. Centring the direction removes it.
    ///
    /// Their worst observed cell was a permuted null running +0.542 to
    /// +0.698 while the real effect was +0.611 — the shuffled arm
    /// separated cases better than the real one.
    #[test]
    fn per_subject_centering_removes_the_anchor_artefact() {
        let train = training();
        // Held-out cohort with NO disease signal anywhere, but whose
        // REFERENCE proteins differ between arms. A uniform bump would
        // not do: reference-protein normalisation subtracts exactly
        // that. The artefact needs the anchor itself to move, which is
        // what happens when the reference set is not arm-neutral in the
        // held-out cohort — and nothing in the fit can know that,
        // because the set was chosen on the training cohorts.
        let mut held = cohort(
            "c",
            Spec {
                n_signal: 0,
                offset: 100.0,
                scale: 7.0,
                ..SPEC
            },
            3,
        );
        let refs: Vec<usize> = match &fit(&train, MethodSpec::ReferenceProtein { k: 8 }, false, 1)
            .unwrap()
            .method
        {
            HarmonizeMethod::ReferenceProtein { reference_features } => reference_features
                .iter()
                .filter_map(|f| held.features.iter().position(|g| g == f))
                .collect(),
            _ => unreachable!(),
        };
        assert!(
            !refs.is_empty(),
            "reference features must map into the held-out cohort"
        );
        for (i, row) in held.values.iter_mut().enumerate() {
            if i % 2 == 0 {
                for j in &refs {
                    row[*j] += 6.0;
                }
            }
        }

        // A direction learned from shuffled labels: no disease content.
        let model = fit(&train, MethodSpec::ReferenceProtein { k: 8 }, true, 5).unwrap();

        let uncentred = separation(&apply(&model, &held).unwrap()).abs();
        let centred =
            separation(&apply_with(&model, &held, DirectionCentering::PerSubject).unwrap()).abs();
        assert!(
            uncentred > 1.0,
            "fixture should reproduce the artefact; got {uncentred:.3}"
        );
        assert!(
            centred < uncentred / 4.0,
            "per-subject centring should remove the anchor term: {uncentred:.3} -> {centred:.3}"
        );
    }

    /// Centring must not destroy a direction that lives in the
    /// contrasts between features rather than in their common level.
    #[test]
    fn per_subject_centering_keeps_a_real_contrast() {
        let train = training();
        let held = held_out();
        let model = fit(&train, MethodSpec::ZScore, false, 42).unwrap();
        let centred =
            separation(&apply_with(&model, &held, DirectionCentering::PerSubject).unwrap()).abs();
        assert!(
            centred > 0.5,
            "a genuine contrast should survive centring; got {centred:.3}"
        );
    }
}
