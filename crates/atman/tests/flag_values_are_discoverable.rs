//! A flag value that the parser accepts but nothing documents is
//! unreachable in practice.
//!
//! This is the least visible member of the defect class this codebase
//! kept hitting. The others produced a wrong or meaningless value that
//! could in principle be caught by reading outputs: an inert threshold,
//! a discarded weight, a truncated `k`. This one produces correct
//! behaviour that cannot be found, with the documentation actively
//! asserting it does not exist. No output can reveal it. It was
//! introduced here by a fix for one of the other instances: the
//! `cosine-centered` metric reached `align bootstrap`'s parser and its
//! rejection message and help text both continued to list three values.
//!
//! Three places must agree for a string-parsed enum flag: what the
//! parser accepts, what the rejection message says it accepts, and what
//! `--help` lists. This test pins all three.

use std::process::{Command, Output};

fn run_atman(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_atman"))
        .args(args)
        .output()
        .expect("run atman")
}

/// `(subcommand path, flag, accepted values)`.
const FLAG_VALUES: &[(&[&str], &str, &[&str])] = &[
    (
        &["align", "bootstrap"],
        "--metric",
        &["cosine", "cosine-centered", "jaccard", "spearman"],
    ),
    (&["align", "bootstrap"], "--decomposition", &["ica", "nmf"]),
    (
        &["decompose", "unmix"],
        "--method",
        &["spa", "vca", "nfindr"],
    ),
    (&["decompose", "unmix"], "--abundance", &["fcls", "ucls"]),
    (
        &["modules", "discover"],
        "--method",
        &["wgcna-soft", "hard-threshold"],
    ),
];

#[test]
fn every_accepted_flag_value_appears_in_help() {
    for (path, flag, values) in FLAG_VALUES {
        let mut args: Vec<&str> = path.to_vec();
        args.push("--help");
        let out = run_atman(&args);
        let help = String::from_utf8_lossy(&out.stdout);
        // Isolate the flag's own help block so a value mentioned under a
        // different flag cannot satisfy the assertion.
        let block = help
            .split(flag)
            .nth(1)
            .unwrap_or_else(|| panic!("{path:?} --help does not mention {flag}"));
        let block = block.split("\n      --").next().unwrap_or(block);
        for v in *values {
            assert!(
                block.contains(v),
                "{path:?} {flag} accepts {v:?} but its help text does not mention it; a value \
                 that works and cannot be found is indistinguishable from one that does not \
                 exist.\nHelp block:\n{block}"
            );
        }
    }
}

/// The rejection message must name every value the parser accepts.
///
/// Checked against the source rather than by invoking the binary: a
/// bogus value cannot reach these `bail!`s from the command line
/// because clap rejects the missing required arguments first, so an
/// invocation-based version of this test skips every case and passes
/// vacuously. That is the same defect it exists to catch, so it is
/// worth stating plainly: the first draft of this test did exactly
/// that, and only a deliberately undocumented canary value revealed it.
#[test]
fn rejection_message_names_every_match_arm() {
    let files = [
        "src/commands/align.rs",
        "src/commands/decompose/unmix.rs",
        "src/commands/decompose/ica.rs",
        "src/commands/modules.rs",
    ];
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut checked = 0usize;
    for rel in files {
        let path = root.join(rel);
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        // Every `other => bail!("--flag ...; expected ...")` and the
        // literal match arms above it in the same match block.
        for (idx, _) in text.match_indices("=> bail!(") {
            let msg_start = text[idx..].find('"').map(|o| idx + o + 1);
            let Some(ms) = msg_start else { continue };
            let Some(me) = text[ms..].find('"').map(|o| ms + o) else {
                continue;
            };
            let message = &text[ms..me];
            if !message.starts_with("--") || !message.contains("expected") {
                continue;
            }
            // Walk back to the opening of this match block, and only
            // accept it if the matched expression names the same flag.
            // Without that check the nearest preceding `match` is often
            // a different flag's, which produced a confident false
            // report on the first run.
            let Some(match_start) = text[..idx].rfind(" match ") else {
                continue;
            };
            let flag_ident = message
                .split_whitespace()
                .next()
                .unwrap_or("")
                .trim_start_matches("--")
                .replace('-', "_");
            let scrutinee_end = text[match_start..idx]
                .find('{')
                .map(|o| match_start + o)
                .unwrap_or(idx);
            if flag_ident.is_empty() || !text[match_start..scrutinee_end].contains(&flag_ident) {
                continue;
            }
            let arms_region = &text[match_start..idx];
            let mut arms: Vec<String> = Vec::new();
            for (ai, _) in arms_region.match_indices("=>") {
                let before = arms_region[..ai].trim_end();
                if !before.ends_with('"') {
                    continue;
                }
                let lit_end = before.len() - 1;
                let Some(lit_start) = before[..lit_end].rfind('"') else {
                    continue;
                };
                let lit = &before[lit_start + 1..lit_end];
                // An arm whose body is itself a `bail!` is a value the
                // parser RECOGNISES in order to refuse it with a
                // specific explanation, which is better than a generic
                // rejection. It is not an accepted value and must not
                // be advertised as one.
                let body = &arms_region[ai + 2..];
                let body_head = &body[..body.len().min(60)];
                let refused = body_head.trim_start().starts_with("bail!")
                    || body_head.trim_start().starts_with("{\n") && body_head.contains("bail!");
                if !lit.is_empty()
                    && !refused
                    && lit.chars().all(|c| c.is_ascii_lowercase() || c == '-')
                {
                    arms.push(lit.to_string());
                }
            }
            if arms.is_empty() {
                continue;
            }
            checked += 1;
            for arm in &arms {
                assert!(
                    message.contains(arm.as_str()),
                    "{rel}: the parser accepts {arm:?} but its rejection message omits it, so a \
                     user told what is valid is told wrong.\n  message: {message}\n  arms: \
                     {arms:?}"
                );
            }
        }
    }
    assert!(
        checked >= 4,
        "expected to check several flags; found {checked}. If the source shape changed this \
         test may be silently checking nothing, which is the failure it exists to prevent."
    );
}
