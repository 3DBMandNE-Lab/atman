use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;
use csv::ReaderBuilder;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::PathBuf;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Per-subject delta table with columns:
    /// subject_id, gene_symbol, pt1_pr1, pt2_pr2, pt2_pt1, pr2_pr1
    #[arg(long)]
    deltas_tsv: PathBuf,

    /// Module definition TSV with columns: module, gene_symbol
    #[arg(long)]
    modules_tsv: PathBuf,

    /// Output path for module trajectory scores TSV.
    #[arg(long)]
    output: PathBuf,
}

#[derive(Default, Clone, Copy)]
struct Acc {
    sum_pt1_pr1: f64,
    sum_pt2_pr2: f64,
    sum_pt2_pt1: f64,
    sum_pr2_pr1: f64,
    n: usize,
}

pub fn run(args: Args) -> Result<()> {
    let modules = read_modules(&args.modules_tsv)?;
    if modules.is_empty() {
        bail!("no module definitions found in {:?}", args.modules_tsv);
    }
    let module_name_by_gene: HashMap<String, Vec<String>> = {
        let mut m: HashMap<String, Vec<String>> = HashMap::new();
        for (module, genes) in &modules {
            for g in genes {
                m.entry(g.clone()).or_default().push(module.clone());
            }
        }
        m
    };

    let mut acc: BTreeMap<(String, String), Acc> = BTreeMap::new();
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&args.deltas_tsv)
        .with_context(|| format!("opening {:?}", args.deltas_tsv))?;
    let headers = rdr.headers()?.clone();
    let need = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing column {:?} in {:?}", name, args.deltas_tsv))
    };
    let c_subject = need("subject_id")?;
    let c_gene = need("gene_symbol")?;
    let c_pt1 = need("pt1_pr1")?;
    let c_pt2 = need("pt2_pr2")?;
    let c_accpre = need("pt2_pt1")?;
    let c_accpost = need("pr2_pr1")?;

    for rec in rdr.records() {
        let row = rec?;
        let gene = row[c_gene].trim();
        let subject = row[c_subject].trim();
        let ms = match module_name_by_gene.get(gene) {
            Some(v) => v,
            None => continue,
        };
        let pt1: f64 = row[c_pt1].trim().parse().with_context(|| "parse pt1_pr1")?;
        let pt2: f64 = row[c_pt2].trim().parse().with_context(|| "parse pt2_pr2")?;
        let ap: f64 = row[c_accpre]
            .trim()
            .parse()
            .with_context(|| "parse pt2_pt1")?;
        let ao: f64 = row[c_accpost]
            .trim()
            .parse()
            .with_context(|| "parse pr2_pr1")?;
        for module in ms {
            let e = acc
                .entry((subject.to_string(), module.clone()))
                .or_default();
            e.sum_pt1_pr1 += pt1;
            e.sum_pt2_pr2 += pt2;
            e.sum_pt2_pt1 += ap;
            e.sum_pr2_pr1 += ao;
            e.n += 1;
        }
    }

    let mut out = String::from("subject_id\tmodule\tn_genes\tpt1_pr1\tpt2_pr2\tpt2_pt1\tpr2_pr1\n");
    for ((subject, module), a) in acc {
        if a.n == 0 {
            continue;
        }
        let n = a.n as f64;
        out.push_str(&subject);
        out.push('\t');
        out.push_str(&module);
        out.push('\t');
        out.push_str(&a.n.to_string());
        out.push('\t');
        out.push_str(&(a.sum_pt1_pr1 / n).to_string());
        out.push('\t');
        out.push_str(&(a.sum_pt2_pr2 / n).to_string());
        out.push('\t');
        out.push_str(&(a.sum_pt2_pt1 / n).to_string());
        out.push('\t');
        out.push_str(&(a.sum_pr2_pr1 / n).to_string());
        out.push('\n');
    }
    std::fs::write(&args.output, out.as_bytes())
        .with_context(|| format!("writing {:?}", args.output))?;
    eprintln!("module-trajectory: wrote {:?}", args.output);
    Ok(())
}

fn read_modules(path: &PathBuf) -> Result<BTreeMap<String, BTreeSet<String>>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)
        .with_context(|| format!("opening {:?}", path))?;
    let headers = rdr.headers()?.clone();
    let need = |name: &str| -> Result<usize> {
        headers
            .iter()
            .position(|h| h == name)
            .with_context(|| format!("missing column {:?} in {:?}", name, path))
    };
    let c_module = need("module")?;
    let c_gene = need("gene_symbol")?;
    let mut out: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for rec in rdr.records() {
        let row = rec?;
        let m = row[c_module].trim();
        let g = row[c_gene].trim();
        if m.is_empty() || g.is_empty() {
            continue;
        }
        out.entry(m.to_string()).or_default().insert(g.to_string());
    }
    Ok(out)
}
