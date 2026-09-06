//! A ratchet on documentation readability.
//!
//! `docs/style.md` states the rules. This measures the two of them a
//! counter can check without judgement: how many sentences run past 25
//! words, and how many semicolons appear. Both correlate with technical
//! prose a reader has to re-read, and both drift upward silently as
//! documents are edited.
//!
//! The budgets record the state on 2026-09-06. They exist to ratchet:
//! new text must not make a file worse. Lower a budget when a file
//! improves, and never raise one.
//!
//! What this cannot see: voice, tense, noun clusters, dropped articles,
//! whether a caveat survived an edit. A file that passes here can still
//! read badly. That is a reason to review prose, not a reason to skip
//! the counter, because review does not catch drift and the counter
//! does.

use std::path::{Path, PathBuf};

/// `(path, max share of long sentences in percent, max semicolons)`.
const BUDGETS: &[(&str, usize, usize)] = &[
    ("README.md", 4, 0),
    ("docs/reference.md", 8, 0),
    ("docs/recipes.md", 17, 14),
    ("docs/tutorial.md", 2, 8),
    ("docs/style.md", 1, 0),
];

const LONG_SENTENCE_WORDS: usize = 25;

fn repo_root() -> PathBuf {
    // CARGO_MANIFEST_DIR is crates/atman.
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

/// Prose only: drop fenced code, table rows, headings, and link-only
/// lines, because none of them is a sentence a reader parses as prose.
fn prose_of(markdown: &str) -> String {
    let mut out = String::new();
    let mut in_fence = false;
    for line in markdown.lines() {
        let t = line.trim_start();
        if t.starts_with("```") {
            in_fence = !in_fence;
            continue;
        }
        if in_fence || t.starts_with('|') || t.starts_with('#') || t.starts_with("[!") {
            continue;
        }
        if t.is_empty() {
            // A blank line ends a paragraph. Without this the extractor
            // glues separate paragraphs and list items into one very
            // long pseudo-sentence, which inflates the long-sentence
            // count and points the failure message at text that is not
            // actually a run-on.
            out.push_str(". ");
            continue;
        }
        out.push_str(line);
        let ends_sentence = t.ends_with('.') || t.ends_with('!') || t.ends_with('?');
        let is_item = t.starts_with("- ")
            || t.starts_with("* ")
            || t.chars().next().is_some_and(|c| c.is_ascii_digit()) && t.contains(". ");
        if is_item && !ends_sentence {
            out.push('.');
        }
        out.push(' ');
    }
    out
}

/// Semicolons that are punctuation, not syntax.
///
/// A semicolon inside a code span or a quoted CLI value is something the
/// reader must type: a `;`-separated flag value, a quoted argument, a
/// MaxQuant protein-group identifier. Counting those would penalise
/// documenting them accurately, which is the opposite of the rule's
/// purpose.
fn prose_semicolons(markdown: &str) -> usize {
    let mut count = 0usize;
    let mut in_fence = false;
    for line in markdown.lines() {
        let t = line.trim_start();
        if t.starts_with("```") {
            in_fence = !in_fence;
            continue;
        }
        if in_fence {
            continue;
        }
        // State resets every line. A markdown code span does not cross a
        // line break, and carrying the state onward means one unbalanced
        // backtick silently changes the count for the whole rest of the
        // file.
        let mut in_code = false;
        let mut in_quote = false;
        for c in line.chars() {
            match c {
                '`' => in_code = !in_code,
                // Double quotes only. An apostrophe in a possessive
                // would otherwise desynchronise the rest of the line.
                '"' => in_quote = !in_quote,
                ';' if !in_code && !in_quote => count += 1,
                _ => {}
            }
        }
    }
    count
}

fn sentences(prose: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut current = String::new();
    let chars: Vec<char> = prose.chars().collect();
    for (i, c) in chars.iter().enumerate() {
        current.push(*c);
        if matches!(c, '.' | '!' | '?') {
            // Not a sentence end inside a version, a decimal, or an
            // abbreviation followed by more of the same word.
            let next_is_space = chars.get(i + 1).map(|n| n.is_whitespace()).unwrap_or(true);
            let prev_is_digit = i > 0 && chars[i - 1].is_ascii_digit();
            let next_is_digit = chars
                .get(i + 1)
                .map(|n| n.is_ascii_digit())
                .unwrap_or(false);
            if next_is_space && !(prev_is_digit && next_is_digit) {
                let s = current.trim().to_string();
                if s.split_whitespace().count() > 2 {
                    out.push(s);
                }
                current.clear();
            }
        }
    }
    let s = current.trim().to_string();
    if s.split_whitespace().count() > 2 {
        out.push(s);
    }
    out
}

#[test]
fn documentation_readability_does_not_regress() {
    let root = repo_root();
    let mut failures: Vec<String> = Vec::new();
    let mut report: Vec<String> = Vec::new();

    for (rel, max_long_pct, max_semicolons) in BUDGETS {
        let path = root.join(rel);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let prose = prose_of(&text);
        let sents = sentences(&prose);
        assert!(
            sents.len() > 5,
            "{rel}: only {} sentences found; the prose extractor is probably broken, which \
             would make every budget below pass for the wrong reason",
            sents.len()
        );
        let long: Vec<&String> = sents
            .iter()
            .filter(|s| s.split_whitespace().count() > LONG_SENTENCE_WORDS)
            .collect();
        let long_pct = 100 * long.len() / sents.len();
        let semicolons = prose_semicolons(&text);
        report.push(format!(
            "  {rel}: {}/{} long ({long_pct}%, budget {max_long_pct}%), {semicolons} semicolons \
             (budget {max_semicolons})",
            long.len(),
            sents.len()
        ));
        if long_pct > *max_long_pct {
            let worst = long
                .iter()
                .max_by_key(|s| s.split_whitespace().count())
                .map(|s| s.split_whitespace().take(18).collect::<Vec<_>>().join(" "))
                .unwrap_or_default();
            failures.push(format!(
                "{rel}: {long_pct}% of sentences exceed {LONG_SENTENCE_WORDS} words, budget is \
                 {max_long_pct}%. Longest starts: \"{worst}...\""
            ));
        }
        if semicolons > *max_semicolons {
            failures.push(format!(
                "{rel}: {semicolons} semicolons, budget is {max_semicolons}. Start a new \
                 sentence instead."
            ));
        }
    }

    if !failures.is_empty() {
        panic!(
            "documentation readability regressed against docs/style.md:\n{}\n\nmeasured:\n{}",
            failures.join("\n"),
            report.join("\n")
        );
    }
}
