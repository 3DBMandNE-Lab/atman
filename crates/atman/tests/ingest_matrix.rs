use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    let bin = env!("CARGO_BIN_EXE_atman");
    Command::new(bin).args(args).output().expect("run atman")
}

#[test]
fn ingest_matrix_supports_proteins_rows_and_preserves_sample_covariates() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tsubject\tgroup\tage\n\
S1\tA\tControl\t30\n\
S2\tB\tCase\t42\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
Protein.Group\tGenes\tPanel\tS1\tS2\n\
P00001\tGENE1\tms\t10\t12\n\
P00002\tGENE2\tms\t8\t\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--orientation",
        "proteins-rows",
        "--platform",
        "spectronaut_report",
        "--abundance-unit",
        "log2_intensity",
        "--sample-id-col",
        "sample_id",
        "--subject-id-col",
        "subject",
        "--condition-col",
        "group",
        "--assay-id-col",
        "Protein.Group",
        "--gene-col",
        "Genes",
        "--panel-col",
        "Panel",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let samples = std::fs::read_to_string(out.join("samples.tsv")).unwrap();
    assert!(samples.starts_with(
        "sample_id\tsubject_id\tcondition\tis_control\tsample_type\tingest_order\tage\n"
    ));
    assert!(samples.contains("S2\tB\tCase\t0\t\t2\t42\n"));

    let proteins = std::fs::read_to_string(out.join("proteins.tsv")).unwrap();
    assert!(proteins.contains("spectronaut_report\tP00001\t\tGENE1\tms\t\n"));

    let measurements = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    assert!(measurements
        .contains("spectronaut_report\tS1\tP00001\tGENE1\tms\t10\t10\t10\tlog2_intensity"));
    assert!(measurements.contains(
        "spectronaut_report\tS2\tP00002\tGENE2\tms\t\t\t\tlog2_intensity\tPASS\tPASS\t\t0\t1"
    ));
    // qc_measurements.tsv is no longer written; the dropped_by_qc flag
    // inside measurements.tsv is the authoritative QC boundary.
    assert!(
        out.join("measurements.tsv").exists(),
        "measurements.tsv must exist after ingest"
    );
    assert!(
        !out.join("qc_measurements.tsv").exists(),
        "qc_measurements.tsv should no longer be written by ingest"
    );
}

#[test]
fn ingest_matrix_supports_samples_rows_and_log2_transform() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tcondition\n\
S1\tControl\n\
S2\tCase\n",
    )
    .unwrap();
    std::fs::write(
        input.join("proteins.tsv"),
        "\
assay\tgene\n\
A1\tGENE1\n\
A2\tGENE2\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
sample_id\tA1\tA2\n\
S1\t8\t0\n\
S2\t16\t4\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--proteins",
        input.join("proteins.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--orientation",
        "samples-rows",
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--gene-col",
        "gene",
        "--log2-transform",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let measurements = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    assert!(measurements.contains("diann_report\tS1\tA1\tGENE1\t\t8\t3\t3\tlog2_diann_pg_quantity"));
    assert!(measurements.contains(
        "diann_report\tS1\tA2\tGENE2\t\t0\t\t\tlog2_diann_pg_quantity\tPASS\tPASS\t\t0\t1"
    ));

    let validate = run_atman(&["validate", "--input-dir", out.to_str().unwrap()]);
    assert!(
        validate.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&validate.stderr)
    );
}

fn parse_measurement_abundance_by_sample(
    contents: &str,
) -> std::collections::HashMap<String, Vec<f64>> {
    let mut lines = contents.lines();
    let header = lines.next().expect("header");
    let cols: Vec<&str> = header.split('\t').collect();
    let sid = cols.iter().position(|c| *c == "sample_id").unwrap();
    let abund = cols.iter().position(|c| *c == "abundance").unwrap();
    let dropped = cols.iter().position(|c| *c == "dropped_by_qc").unwrap();
    let mut out: std::collections::HashMap<String, Vec<f64>> = Default::default();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.get(dropped).copied() == Some("1") {
            continue;
        }
        let Some(v) = fields
            .get(abund)
            .and_then(|s| (!s.is_empty()).then(|| s.parse::<f64>().ok()).flatten())
        else {
            continue;
        };
        out.entry(fields[sid].to_string()).or_default().push(v);
    }
    out
}

fn sample_median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    }
}

/// Parse one named numeric column of measurements.tsv, grouped by sample_id,
/// skipping QC-dropped rows. Used to compare `abundance` vs `abundance_raw`.
fn parse_named_column_by_sample(
    contents: &str,
    col_name: &str,
) -> std::collections::HashMap<String, Vec<f64>> {
    let mut lines = contents.lines();
    let header = lines.next().expect("header");
    let cols: Vec<&str> = header.split('\t').collect();
    let sid = cols.iter().position(|c| *c == "sample_id").unwrap();
    let col = cols
        .iter()
        .position(|c| *c == col_name)
        .unwrap_or_else(|| panic!("missing column {col_name}"));
    let dropped = cols.iter().position(|c| *c == "dropped_by_qc").unwrap();
    let mut out: std::collections::HashMap<String, Vec<f64>> = Default::default();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.get(dropped).copied() == Some("1") {
            continue;
        }
        let Some(v) = fields
            .get(col)
            .and_then(|s| (!s.is_empty()).then(|| s.parse::<f64>().ok()).flatten())
        else {
            continue;
        };
        out.entry(fields[sid].to_string()).or_default().push(v);
    }
    out
}

/// TASK-001 fidelity: `abundance_raw` must hold the pre-normalization input
/// value, NOT the normalized `abundance`. With S1 shifted up 4 log2 units,
/// median normalization moves `abundance` but must leave `abundance_raw`
/// equal to the raw inputs {14,16,18,20}.
#[test]
fn ingest_matrix_abundance_raw_is_pre_normalization() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tcondition\nS1\tControl\nS2\tCase\nS3\tControl\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tS1\tS2\tS3\n\
A1\t14\t10\t12\n\
A2\t16\t12\t14\n\
A3\t18\t14\t16\n\
A4\t20\t16\t18\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--normalize",
        "median",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let contents = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    let raw = parse_named_column_by_sample(&contents, "abundance_raw");
    let norm = parse_named_column_by_sample(&contents, "abundance");

    // abundance_raw for S1 is exactly the raw inputs, untouched by normalization.
    let mut s1_raw = raw.get("S1").cloned().unwrap();
    s1_raw.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(
        s1_raw,
        vec![14.0, 16.0, 18.0, 20.0],
        "abundance_raw mutated"
    );

    // Median normalization shifts S1's working abundance, so the two columns
    // must differ for at least one row — proving abundance_raw is not a copy
    // of the normalized value (the original bug).
    let s1_norm = norm.get("S1").unwrap();
    let s1_raw_unsorted = raw.get("S1").unwrap();
    assert!(
        s1_norm
            .iter()
            .zip(s1_raw_unsorted)
            .any(|(n, r)| (n - r).abs() > 1e-9),
        "abundance and abundance_raw identical after normalization — raw not preserved"
    );
}

#[test]
fn ingest_matrix_median_normalization_equalizes_sample_medians() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tcondition\nS1\tControl\nS2\tCase\nS3\tControl\n",
    )
    .unwrap();
    // Sample S1 is shifted up by 4 log2 units vs S2; S3 sits in between.
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tS1\tS2\tS3\n\
A1\t14\t10\t12\n\
A2\t16\t12\t14\n\
A3\t18\t14\t16\n\
A4\t20\t16\t18\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--normalize",
        "median",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let contents = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    let mut by_sample = parse_measurement_abundance_by_sample(&contents);
    let s1 = sample_median(by_sample.get_mut("S1").unwrap());
    let s2 = sample_median(by_sample.get_mut("S2").unwrap());
    let s3 = sample_median(by_sample.get_mut("S3").unwrap());
    assert!((s1 - s2).abs() < 1e-9, "S1 {s1} S2 {s2}");
    assert!((s1 - s3).abs() < 1e-9, "S1 {s1} S3 {s3}");
}

#[test]
fn ingest_matrix_quantile_normalization_matches_distributions() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tcondition\nS1\tControl\nS2\tCase\n",
    )
    .unwrap();
    // S1 and S2 have different value ranges and different orderings.
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tS1\tS2\n\
A1\t1\t10\n\
A2\t2\t8\n\
A3\t3\t6\n\
A4\t4\t4\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--normalize",
        "quantile",
    ]);
    assert!(
        output.status.success(),
        "stderr:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let contents = std::fs::read_to_string(out.join("measurements.tsv")).unwrap();
    let mut by_sample = parse_measurement_abundance_by_sample(&contents);
    let mut s1 = by_sample.remove("S1").unwrap();
    let mut s2 = by_sample.remove("S2").unwrap();
    s1.sort_by(|a, b| a.partial_cmp(b).unwrap());
    s2.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(s1.len(), s2.len());
    for (a, b) in s1.iter().zip(s2.iter()) {
        assert!(
            (a - b).abs() < 1e-9,
            "sorted distributions differ: {a} vs {b}"
        );
    }
}

#[test]
fn ingest_matrix_warns_when_sample_medians_diverge_without_normalization() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "sample_id\tcondition\nS1\tControl\nS2\tCase\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tS1\tS2\n\
A1\t14\t10\n\
A2\t16\t12\n\
A3\t18\t14\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "diann_report",
        "--abundance-unit",
        "log2_diann_pg_quantity",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
    ]);
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("per-sample median abundance spans"),
        "expected median-spread warning; stderr:\n{stderr}"
    );
}

#[test]
fn ingest_matrix_fails_when_no_sample_columns_match() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let out = input.join("out");
    std::fs::write(
        input.join("samples.tsv"),
        "\
sample_id\tcondition\n\
S1\tControl\n",
    )
    .unwrap();
    std::fs::write(
        input.join("matrix.tsv"),
        "\
assay\tgene\tX1\n\
A1\tGENE1\t1\n",
    )
    .unwrap();

    let output = run_atman(&[
        "ingest-matrix",
        "--matrix",
        input.join("matrix.tsv").to_str().unwrap(),
        "--samples",
        input.join("samples.tsv").to_str().unwrap(),
        "--output-dir",
        out.to_str().unwrap(),
        "--platform",
        "maxquant_lfq",
        "--abundance-unit",
        "log2_lfq",
        "--condition-col",
        "condition",
        "--assay-id-col",
        "assay",
        "--gene-col",
        "gene",
    ]);

    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("no matrix columns matched sample IDs")
    );
}

/// `--normalize median` must be bit-reproducible across invocations of the
/// same binary on the same input. Every cell of `abundance` carries the grand
/// mean of the per-sample medians, so the order in which those medians are
/// summed must not depend on anything that varies between processes.
///
/// Sixty-four samples with log2-transformed linear intensities give medians
/// whose low bits differ, so a summation order that varies between runs
/// changes the grand mean in the last place and shows up as a different file.
#[test]
fn ingest_matrix_median_normalization_is_bit_reproducible_across_invocations() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path();
    let n_samples = 64;
    let n_assays = 40;

    let mut samples = String::from("sample_id\tcondition\n");
    for s in 0..n_samples {
        samples.push_str(&format!(
            "S{s:03}\t{}\n",
            if s % 2 == 0 { "Control" } else { "Case" }
        ));
    }
    std::fs::write(input.join("samples.tsv"), samples).unwrap();

    // Deterministic pseudo-random linear intensities (LCG); no external RNG.
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let unit = ((state >> 11) as f64) / ((1u64 << 53) as f64);
        1000.0 + unit * 2.0e6
    };
    let mut matrix = String::from("assay");
    for s in 0..n_samples {
        matrix.push_str(&format!("\tS{s:03}"));
    }
    matrix.push('\n');
    for a in 0..n_assays {
        matrix.push_str(&format!("P{a:05}"));
        for _ in 0..n_samples {
            matrix.push_str(&format!("\t{:.6}", next()));
        }
        matrix.push('\n');
    }
    std::fs::write(input.join("matrix.tsv"), matrix).unwrap();

    let mut outputs: Vec<Vec<u8>> = Vec::new();
    for run in 0..6 {
        let out = input.join(format!("run{run}"));
        let output = run_atman(&[
            "ingest-matrix",
            "--matrix",
            input.join("matrix.tsv").to_str().unwrap(),
            "--samples",
            input.join("samples.tsv").to_str().unwrap(),
            "--output-dir",
            out.to_str().unwrap(),
            "--platform",
            "diann_report",
            "--abundance-unit",
            "linear_pg_quantity",
            "--condition-col",
            "condition",
            "--assay-id-col",
            "assay",
            "--log2-transform",
            "--normalize",
            "median",
        ]);
        assert!(
            output.status.success(),
            "run {run} stderr:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        outputs.push(std::fs::read(out.join("measurements.tsv")).unwrap());
    }

    for (run, bytes) in outputs.iter().enumerate().skip(1) {
        assert!(
            bytes == &outputs[0],
            "measurements.tsv from run {run} differs from run 0: \
median normalization is not reproducible across invocations"
        );
    }
}
