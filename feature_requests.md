# Atman feature requests

Inbox for feature proposals **a project has actually asked for**.
Shipped work is documented in `CHANGELOG.md`; this file tracks only
what has been requested and is not yet shipped.

An entry belongs here only if a named requester needs it for real work.
Speculative method-transfer ideas go to the "Unscheduled" section of
`docs/analytical-roadmap.md` instead — an inbox where most entries are
nobody's request teaches the next session to ignore it.

Add new proposals as sections below, naming the requester. Remove
entries once they land or once the need goes away (the commit message
and `CHANGELOG.md` are the permanent record).

---

## Open

### Hierarchical / nested alignment

Site → cohort → cross-cohort alignment for consortium data where a
single cohort has multiple acquisition sites. Each level of the
hierarchy carries its own bootstrap-based confidence interval
(building on the existing `atman align bootstrap` surface).

**Designed, not scheduled.** Spec at
`docs/superpowers/specs/2026-09-06-hierarchical-alignment-design.md`.
Requested by the CSF CrossDisease session as a revision tool (a
reviewer asking whether Axis 1 replicates within Bader's four
collection sites); not needed for that paper's submission, and no
other manuscript session needs it. Build it when a reviewer asks.

---

## Demand survey, 2026-09-06

Every manuscript session on this machine was asked which open item it
needs. Meningioma (Text/v4) and CAR-T do not use atman at all. The GBM
methods paper needs nothing: its run tree is complete at PIN2 and
replays from a snapshotted binary. The CSF CrossDisease paper needs
nothing for submission.

That survey emptied most of this file. Every remaining entry from the
2026-04-20 batch was a speculative port from an adjacent field rather
than anyone's request, so:

- **Pooled-QC drift correction — removed, not deferred.** It has no
  possible input. No cohort in any active project records injection
  order, acquisition date, or pooled-QC samples, and the SIH `plate_id`
  and `panel_lot` columns are empty. `ingest_order` is row order at
  ingest, not acquisition order. Revive it only if data that carries
  acquisition order ever arrives.
- **Compositional effect size for DE — removed as contradicted.** A
  `--effect-size-scale clr` flag conflicts with the stated Methods of
  the only paper it would serve (log2 + median normalization, Cohen d on
  that scale). It would describe a different pipeline, not an option
  within this one.
- **PMF, MCR-ALS and longitudinal tensor decomposition — moved** to the
  "Unscheduled: alternative decompositions" section of
  `docs/analytical-roadmap.md`. Each competes with a decomposition
  approach the current projects have committed to, so each would be a
  new sensitivity analysis nobody asked for. They are ideas, not
  requests, and the roadmap is where ideas live.

Hierarchical alignment stays because a real question sits behind it and
it now has a spec. Re-ask before building anything.
