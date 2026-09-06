//! Module-level differential abundance. Aggregates proteins into
//! data-driven module scores (mean abundance across member proteins per
//! sample) and tests at module level. Moves the testing burden from P
//! proteins to K modules, reducing the BH correction denominator by orders
//! of magnitude.

use anyhow::{bail, Context, Result};
use atman_core::de::{bh_fdr, paired_t, welch_t, PairedTResult, SkipReason};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::time::SystemTime;

use super::parse_comparisons;
use crate::io::{
    hash_canonical_inputs, hash_labeled_inputs, read_measurements_long, sidecar_path_for,
    write_run_sidecar,
};

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Directory containing measurements.tsv and samples.tsv.
    #[arg(long)]
    input_dir: PathBuf,

    /// Output directory for module_de_results.tsv.
    #[arg(long)]
    output_dir: PathBuf,

    /// Modules TSV with columns `module` and `gene_symbol`.
    #[arg(long)]
    modules_tsv: PathBuf,

    /// Test kind: `paired-t` (default) or `welch-t`.
    #[arg(long, default_value = "paired-t")]
    test: String,

    /// Biological-replicate key. Currently fixed to the canonical
    /// `subject_id` column from samples.tsv. Ignored for welch-t.
    #[arg(long, default_value = "subject_id")]
    paired_by: String,

    /// Comma-separated comparisons in `A-B` form.
    #[arg(long)]
    groups: String,

    /// Minimum subjects per group for a test to run.
    #[arg(long, default_value_t = 5)]
    min_pairs: usize,
}

pub fn run(args: Args) -> Result<()> {
    let started_at = SystemTime::now();
    if args.test != "paired-t" && args.test != "welch-t" {
        bail!(
            "test {:?} not supported here; use paired-t or welch-t",
            args.test
        );
    }
    let is_unpaired = args.test == "welch-t";
    if !is_unpaired && args.paired_by != "subject_id" {
        bail!(
            "paired-by {:?} not supported; pairing currently uses the canonical \
             `subject_id` column from samples.tsv.",
            args.paired_by
        );
    }
    std::fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {:?}", args.output_dir))?;

    let comparisons = parse_comparisons(&args.groups)?;
    let modules = read_modules(&args.modules_tsv)?;
    // A module set whose only member is the grey catch-all is not a
    // module set: every gene is in one bucket, so the "per-module" test
    // is a single test of the mean of every feature, which is
    // well-formed, meaningless, and can easily come out significant.
    // `modules discover` refuses to emit such a set, but this file may
    // predate that check or come from elsewhere.
    let non_grey: Vec<&String> = modules.keys().filter(|m| m.as_str() != "grey").collect();
    if non_grey.is_empty() {
        anyhow::bail!(
            "module-de: {:?} contains no module other than the grey catch-all ({} genes). A \
             single all-features module is not a module: the resulting test compares the mean \
             of every feature between groups and its p-value means nothing. Re-run `modules \
             discover` with a lower --cut-height, a smaller --min-module-size, or an explicit \
             --soft-power.",
            args.modules_tsv,
            modules.get("grey").map(|g| g.len()).unwrap_or(0)
        );
    }
    let measurements = read_measurements_long(&args.input_dir.join("measurements.tsv"))?;
    let samples = crate::io::read_samples(&args.input_dir.join("samples.tsv"))?;

    let sample_by_id: HashMap<&str, &atman_core::Sample> =
        samples.iter().map(|s| (s.sample_id.as_str(), s)).collect();

    // gene → list of modules (a gene typically belongs to one module; if not,
    // it contributes to each it's assigned to).
    let mut gene_to_modules: HashMap<&str, Vec<&str>> = HashMap::new();
    for (module, genes) in &modules {
        for g in genes {
            gene_to_modules
                .entry(g.as_str())
                .or_default()
                .push(module.as_str());
        }
    }

    // Accumulate per-(module, sample) abundance sum + count across contributing proteins.
    // Only rows with effective_abundance() (i.e., not QC-masked, not below_lod) contribute.
    #[derive(Default)]
    struct Acc {
        sum: f64,
        n: usize,
    }
    let mut per_module_sample: HashMap<(String, String), Acc> = HashMap::new();
    for m in &measurements {
        let gene = match &m.gene_symbol {
            Some(g) => g.as_str(),
            None => continue,
        };
        let module_list = match gene_to_modules.get(gene) {
            Some(v) => v,
            None => continue,
        };
        let abundance = match m.effective_abundance() {
            Some(a) => a,
            None => continue,
        };
        let sid = m.sample_id.as_str();
        for &mod_name in module_list {
            let entry = per_module_sample
                .entry((mod_name.to_string(), sid.to_string()))
                .or_default();
            entry.sum += abundance;
            entry.n += 1;
        }
    }

    // Build per-(module, condition) → Vec<(subject, module_score)> similar to de.rs.
    type CellsByCondition = BTreeMap<String, Vec<(String, f64)>>;
    let mut cells: BTreeMap<String, CellsByCondition> = BTreeMap::new();
    for ((module, sample_id), acc) in &per_module_sample {
        if acc.n == 0 {
            continue;
        }
        let s = match sample_by_id.get(sample_id.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        if s.is_control {
            continue;
        }
        let subject = match &s.subject_id {
            Some(id) => id.clone(),
            None => continue,
        };
        let condition = match &s.condition {
            Some(c) => c.clone(),
            None => continue,
        };
        let module_score = acc.sum / (acc.n as f64);
        cells
            .entry(module.clone())
            .or_default()
            .entry(condition)
            .or_default()
            .push((subject, module_score));
    }

    // n_genes per module (for the output row; count genes that appear in the measurements).
    let mut module_n_genes: BTreeMap<String, usize> = BTreeMap::new();
    for (module, genes) in &modules {
        let n = genes
            .iter()
            .filter(|g| {
                measurements
                    .iter()
                    .any(|m| m.gene_symbol.as_deref() == Some(g.as_str()))
            })
            .count();
        module_n_genes.insert(module.clone(), n);
    }

    let mut out = String::from(
        "module\tcontrast\tn_samples_a\tn_samples_b\tn_genes\tmean_a\tmean_b\tmean_diff\tt\tdf\tp_value\tbh_q\tskip_reason\n",
    );

    for (comp_a, comp_b) in &comparisons {
        let comparison_label = format!("{}-{}", comp_a, comp_b);
        let mut family_p: Vec<Option<f64>> = Vec::new();
        let mut family_rows: Vec<(String, String, usize, usize, usize, PairedTResult)> = Vec::new();

        for (module, by_cond) in &cells {
            let empty: Vec<(String, f64)> = Vec::new();
            let va = by_cond.get(comp_a).unwrap_or(&empty);
            let vb = by_cond.get(comp_b).unwrap_or(&empty);
            let na = va.len();
            let nb = vb.len();
            let n_genes = *module_n_genes.get(module).unwrap_or(&0);

            let result = if is_unpaired {
                let a_vals: Vec<f64> = va.iter().map(|(_, v)| *v).collect();
                let b_vals: Vec<f64> = vb.iter().map(|(_, v)| *v).collect();
                welch_t(&a_vals, &b_vals, args.min_pairs)
            } else {
                let a_by_subj: HashMap<&str, f64> =
                    va.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                let b_by_subj: HashMap<&str, f64> =
                    vb.iter().map(|(s, v)| (s.as_str(), *v)).collect();
                let pairs: Vec<(f64, f64)> = a_by_subj
                    .iter()
                    .filter_map(|(subj, a)| b_by_subj.get(subj).map(|b| (*a, *b)))
                    .collect();
                paired_t(&pairs, args.min_pairs)
            };
            let p = match &result {
                PairedTResult::Computed { p_value, .. } => Some(*p_value),
                _ => None,
            };
            family_p.push(p);
            family_rows.push((
                module.clone(),
                comparison_label.clone(),
                na,
                nb,
                n_genes,
                result,
            ));
        }
        let qs = bh_fdr(&family_p);
        for ((module, lbl, na, nb, n_genes, r), q) in family_rows.into_iter().zip(qs.into_iter()) {
            match r {
                PairedTResult::Computed {
                    mean_a,
                    mean_b,
                    mean_diff,
                    t,
                    df,
                    p_value,
                    ..
                } => {
                    out.push_str(&format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\n",
                        module,
                        lbl,
                        na,
                        nb,
                        n_genes,
                        mean_a,
                        mean_b,
                        mean_diff,
                        t,
                        df,
                        p_value,
                        q.map(|x| x.to_string()).unwrap_or_default(),
                    ));
                }
                PairedTResult::Skipped { reason, .. } => {
                    let reason_str = match reason {
                        SkipReason::InsufficientPairs => "insufficient_pairs",
                        SkipReason::ZeroVariance => "zero_variance",
                        SkipReason::NonFiniteInput => "non_finite_input",
                    };
                    out.push_str(&format!(
                        "{}\t{}\t{}\t{}\t{}\t\t\t\t\t\t\t\t{}\n",
                        module, lbl, na, nb, n_genes, reason_str
                    ));
                }
            }
        }
    }

    let output_path = args.output_dir.join("module_de_results.tsv");
    crate::io::atomic_write(&output_path, out.as_bytes())
        .with_context(|| "writing module_de_results.tsv")?;
    eprintln!(
        "module-de: wrote {:?} (K={} modules × {} comparisons)",
        output_path,
        cells.len(),
        comparisons.len(),
    );

    let finished_at = SystemTime::now();
    let mut canonical = hash_canonical_inputs(
        &args.input_dir,
        &["measurements.tsv", "samples.tsv", "proteins.tsv"],
    )?;
    let modules_hash = hash_labeled_inputs(&[("modules_tsv", args.modules_tsv.as_path())])?;
    canonical.extend(modules_hash);
    let sidecar = sidecar_path_for(&output_path);
    write_run_sidecar(
        &sidecar,
        "module-de",
        json!({
            "input-dir": args.input_dir.display().to_string(),
            "output-dir": args.output_dir.display().to_string(),
            "modules-tsv": args.modules_tsv.display().to_string(),
            "test": args.test,
            "paired-by": args.paired_by,
            "groups": args.groups,
            "min-pairs": args.min_pairs,
        }),
        &canonical,
        std::slice::from_ref(&output_path),
        started_at,
        finished_at,
        None,
    )?;
    eprintln!("module-de: sidecar={}", sidecar.display());
    Ok(())
}

fn read_modules(path: &PathBuf) -> Result<BTreeMap<String, Vec<String>>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let i_module = headers
        .iter()
        .position(|h| h == "module")
        .with_context(|| format!("missing column `module` in {:?}", path))?;
    let i_gene = headers
        .iter()
        .position(|h| h == "gene_symbol")
        .with_context(|| format!("missing column `gene_symbol` in {:?}", path))?;
    let mut out: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for rec in rdr.records() {
        let row = rec?;
        let m = row[i_module].trim().to_string();
        let g = row[i_gene].trim().to_string();
        if m.is_empty() || g.is_empty() {
            continue;
        }
        out.entry(m).or_default().push(g);
    }
    Ok(out)
}
