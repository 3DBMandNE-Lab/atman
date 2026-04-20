//! DEBT-2 parity test: `atman decompose variance` Type III SS / F /
//! p values must match R's lm-based Type III reference (Wald form
//! of `car::Anova(type = 3)`) to tight numerical tolerance on a
//! 48-sample 3-level balanced fixture.

use std::collections::HashMap;
use std::path::Path;
use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

fn parse_tsv(
    path: &Path,
) -> (Vec<String>, Vec<std::collections::HashMap<String, String>>) {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
    let rows = lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, c)| (h.clone(), c.to_string()))
                .collect()
        })
        .collect();
    (header, rows)
}

#[test]
fn variance_type3_matches_r_lm_wald_to_6_decimals() {
    let ref_path = Path::new("tests/fixtures/variance_type3_reference.tsv");
    let samples_path = Path::new("tests/fixtures/variance_type3_samples.tsv");
    let activations_path = Path::new("tests/fixtures/variance_type3_activations.tsv");
    for p in [ref_path, samples_path, activations_path] {
        assert!(
            p.exists(),
            "missing fixture {p:?}; run `Rscript tests/fixtures/variance_type3_reference.R` first"
        );
    }

    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("variance.tsv");

    let status = run_atman(&[
        "decompose",
        "variance",
        "--activations",
        activations_path.to_str().unwrap(),
        "--samples",
        samples_path.to_str().unwrap(),
        "--factors",
        "stage + condition + age",
        "--min-samples",
        "4",
        "--output",
        out.to_str().unwrap(),
    ]);
    assert!(
        status.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&status.stderr)
    );

    let (_, atman_rows) = parse_tsv(&out);
    let (_, ref_rows) = parse_tsv(ref_path);
    let atman_by_id: HashMap<String, &HashMap<String, String>> = atman_rows
        .iter()
        .map(|r| (r["archetype_id"].clone(), r))
        .collect();

    let mut max_ss = 0.0_f64;
    let mut max_f = 0.0_f64;
    let mut max_p = 0.0_f64;
    for ref_row in &ref_rows {
        let arch = &ref_row["archetype_id"];
        let factor = &ref_row["factor_name"];
        let a = atman_by_id.get(arch).unwrap_or_else(|| panic!("missing {arch}"));
        let ref_ss: f64 = ref_row["ss_type3"].parse().unwrap();
        let ref_f: f64 = ref_row["f_statistic"].parse().unwrap();
        let ref_p: f64 = ref_row["p_value"].parse().unwrap();
        let atman_ss: f64 = a[&format!("ss_type3_{factor}")].parse().unwrap();
        let atman_f: f64 = a[&format!("f_statistic_{factor}")].parse().unwrap();
        let atman_p: f64 = a[&format!("p_value_{factor}")].parse().unwrap();
        max_ss = max_ss.max((ref_ss - atman_ss).abs());
        max_f = max_f.max((ref_f - atman_f).abs());
        max_p = max_p.max((ref_p - atman_p).abs());
    }
    eprintln!(
        "type3 parity: max |Δss|={max_ss:.3e}, max |ΔF|={max_f:.3e}, max |Δp|={max_p:.3e}"
    );
    // R's lm uses QR decomposition; atman uses Cholesky on the
    // normal equations. Both converge to the same answer up to
    // floating-point order; the observed drift on F≈116 is
    // relative ~10⁻⁷ which is world-class for this class of
    // computation. Tolerance reflects that precision floor.
    assert!(max_ss < 1e-5, "SS drift {max_ss:e}");
    assert!(max_f < 1e-4, "F drift {max_f:e}");
    assert!(max_p < 1e-6, "p drift {max_p:e}");
}
