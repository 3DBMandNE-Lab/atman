# Documentation style

Atman's documentation follows ASD-STE100 Simplified Technical English.
STE is a controlled English for aircraft maintenance manuals. It exists
because a mechanic must read a procedure once, under time pressure, and
act on it correctly. That is close to the situation of a user who runs a
command against a cohort that took three years to collect.

Two modes apply here. Pick the mode from what the text does, not from
which file it is in.

## Strict mode

Use strict mode for text a reader acts on directly:

- command-line help (`#[arg(long, ...)]` doc comments)
- error and warning messages
- procedures and step lists in `docs/recipes.md`
- the flag tables in `docs/reference.md`

Rules:

1. One instruction per sentence.
2. Maximum 20 words in a sentence.
3. Active voice. Name the actor. Write "the command refuses" and not
   "an error is raised".
4. Present tense. Do not write present perfect.
5. No semicolons. Start a new sentence.
6. One word, one meaning. A flag is always a "flag" and never also an
   "option", a "switch", or a "parameter".
7. No phrasal verbs where one word does the job. Write "remove" and not
   "take out".
8. No nominalised verbs. Write "the run fails" and not "failure of the
   run occurs".
9. Maximum three nouns in a row. Break up "run sidecar hash column".
10. Keep articles. Write "the sidecar records the seed" and not
    "sidecar records seed".
11. State the condition before the action. Write "if the fit does not
    converge, the command warns" and not the reverse.
12. Keep every fact, condition, and scope qualifier. Brevity never
    removes a caveat.

## STE-flavoured mode

Use STE-flavoured mode for text a reader thinks with:

- `README.md`
- `docs/tutorial.md` prose
- rationale paragraphs in `docs/reference.md`
- `CHANGELOG.md` entries

Keep rules 1 to 5 and rule 12. Relax the vocabulary rules. A rationale
paragraph may name a statistical concept precisely even when the word is
not plain, because the alternative is imprecision.

## What this is not

STE does not mean short documentation. Rule 12 outranks every other
rule. A long explanation written in short, active, single-clause
sentences is correct STE. A brief sentence that drops a condition is
not.

STE also does not mean simple content. The reader of these documents
knows proteomics. They may not know which of two flags controls the
thing that just went wrong.

## Measured conformance

`crates/atman/tests/docs_style.rs` measures two properties that a script
can count without judgement: the share of sentences longer than 25
words, and the number of semicolons. It asserts a per-file budget.

The budgets record the state on 2026-09-07 and exist to ratchet. A new
document must not make a file worse. Lower a budget when you improve a
file, and never raise one.

| File | Long sentences | Semicolons |
|---|---|---|
| `README.md` | 2% | 0 |
| `docs/reference.md` | 7% | 0 |
| `docs/recipes.md` | 17% | 14 |
| `docs/tutorial.md` | 2% | 8 |
| `docs/style.md` | 0% | 0 |

The counter ignores a semicolon inside a code span or a quoted value. A
`;`-separated flag value, a quoted CLI argument, and a MaxQuant
protein-group identifier are all syntax the reader must type. Counting
them would penalise documenting them accurately.

The reference document and the README now carry no prose semicolon. The
recipes and the tutorial still do, and they are the next target.

The counter is deliberately crude. It strips fenced code and table rows,
then splits on sentence-ending punctuation. It cannot see voice, tense,
or noun clusters. A file that passes the counter can still be bad STE.
The counter catches the two failures that correlate best with unreadable
technical prose, and it catches drift, which review does not.
